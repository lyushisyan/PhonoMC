#include "solver/RtaCollisionOperator.h"

#include "PhononMaterial.h"
#include "solver/ParticleStorage.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace phonomc {
namespace {
double lifetime_temperature(double local, const CollisionStep& step) {
    if (!step.local_lifetime_temperature) return step.lifetime_temperature_kelvin;
    return std::isfinite(local) && local > 0.0 ? local : 300.0;
}
}  // namespace

double RtaCollisionOperator::apply(
    ParticleStorage& particles,
    const std::vector<const PhononMaterial*>& materials,
    const CollisionGridView& grid,
    const CollisionStep& step) {
    const int nsv = static_cast<int>(grid.particle_counts.size());
    if (particles.size() <= 0 || nsv <= 0) {
        return 0.0;
    }

    assert(particles.aligned());
    assert(grid.material_ids.size() == grid.particle_counts.size());
    assert(grid.temperatures.size() == grid.particle_counts.size());
    assert(std::accumulate(grid.particle_counts.begin(), grid.particle_counts.end(), 0) == particles.size());

    const CellParticleIndex* cells = grid.cells;
    if (cells == nullptr) {
        fallback_cells_.invalidate();
        fallback_cells_.ensure(particles, nsv);
        cells = &fallback_cells_;
    }
    if (!cells->matches(particles, nsv) || cells->counts() != grid.particle_counts)
        throw std::invalid_argument("RTA requires a current cell index and matching counts.");
    const auto& offsets = cells->offsets();
    const auto& particle_indices = cells->indices();

    auto& relaxation_weight = weights_;
    relaxation_weight.resize(static_cast<size_t>(particles.size()));
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particles.size(); ++i) {
        const PhononMaterial& phonon = *materials[static_cast<size_t>(particles.material_ids[i])];
        relaxation_weight[static_cast<size_t>(i)] = phonon.mode_relaxation_weight(
            lifetime_temperature(particles.temperatures[i], step), particles.modes[i], step.time_step_ps);
    }

    constexpr double kHbarEvPs = 6.582119569e-4;
    constexpr double kBoltzmannEvK = 8.617333262145e-5;
    double residual_raw_ev = 0.0;
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for reduction(+:residual_raw_ev) schedule(dynamic, 16)
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const PhononMaterial& phonon = *materials[static_cast<size_t>(grid.material_ids[static_cast<size_t>(sv)])];
        const double cell_energy_weight = static_cast<double>(phonon.active_mode_count()) *
            particles.spatial_volume(grid.material_ids[sv],step.particle_volume_a3) /
            phonon.energy_density_normalization();
        const int begin = offsets[static_cast<size_t>(sv)];
        const int end = offsets[static_cast<size_t>(sv + 1)];
        if (end <= begin) {
            continue;
        }
        auto evaluate = [&](double temperature, double* derivative) {
            double value = 0.0;
            double slope = 0.0;
            const double T = std::max(0.0, temperature);
            for (int pos = begin; pos < end; ++pos) {
                const int i = particle_indices[static_cast<size_t>(pos)];
                const double omega = std::max(0.0, phonon.mode_angular_frequency(particles.modes[static_cast<size_t>(i)]));
                const double hw = kHbarEvPs * omega;
                const double a = relaxation_weight[static_cast<size_t>(i)];
                if (!(hw > 0.0) || !(a > 0.0)) {
                    continue;
                }
                const double neq = phonon.bose_occupation(T, particles.modes[static_cast<size_t>(i)]);
                value += hw * a * (particles.occupation[static_cast<size_t>(i)] - neq);
                if (T > 0.0) {
                    const double dndt = hw * neq * (neq + 1.0) / (kBoltzmannEvK * T * T);
                    slope -= hw * a * dndt;
                }
            }
            if (derivative != nullptr) {
                *derivative = slope;
            }
            return value;
        };

        double low = 0.0;
        double high = std::max(1000.0, 2.0 * std::max(1.0, grid.temperatures[static_cast<size_t>(sv)]));
        const double f_low = evaluate(low, nullptr);
        double f_high = evaluate(high, nullptr);
        for (int expand = 0; f_high > 0.0 && expand < 20; ++expand) {
            high *= 2.0;
            f_high = evaluate(high, nullptr);
        }

        double pseudo_temperature = std::clamp(grid.temperatures[static_cast<size_t>(sv)], low, high);
        if (f_low <= 0.0) {
            pseudo_temperature = low;
        } else if (f_high >= 0.0) {
            pseudo_temperature = high;
        } else {
            double energy_scale = 0.0;
            for (int pos = begin; pos < end; ++pos) {
                const int i = particle_indices[static_cast<size_t>(pos)];
                const double hw = kHbarEvPs * std::max(
                    0.0, phonon.mode_angular_frequency(particles.modes[static_cast<size_t>(i)]));
                energy_scale += hw * relaxation_weight[static_cast<size_t>(i)] *
                    std::max(1.0, particles.occupation[static_cast<size_t>(i)]);
            }
            for (int iter = 0; iter < 40; ++iter) {
                double derivative = 0.0;
                const double value = evaluate(pseudo_temperature, &derivative);
                if (std::abs(value) <= 1e-13 * std::max(1.0, energy_scale)) {
                    break;
                }
                if (value > 0.0) {
                    low = pseudo_temperature;
                } else {
                    high = pseudo_temperature;
                }
                double candidate = 0.5 * (low + high);
                if (derivative < 0.0 && std::isfinite(derivative)) {
                    const double newton = pseudo_temperature - value / derivative;
                    if (newton > low && newton < high && std::isfinite(newton)) {
                        candidate = newton;
                    }
                }
                pseudo_temperature = candidate;
            }
        }

        double energy_before = 0.0;
        double energy_after = 0.0;
        for (int pos = begin; pos < end; ++pos) {
            const int i = particle_indices[static_cast<size_t>(pos)];
            const double omega = std::max(0.0, phonon.mode_angular_frequency(particles.modes[static_cast<size_t>(i)]));
            const double hw = kHbarEvPs * omega;
            const double a = relaxation_weight[static_cast<size_t>(i)];
            const double old_occupation = particles.occupation[static_cast<size_t>(i)];
            const double neq = phonon.bose_occupation(pseudo_temperature, particles.modes[static_cast<size_t>(i)]);
            const double new_occupation = old_occupation + a * (neq - old_occupation);
            energy_before += hw * old_occupation;
            energy_after += hw * new_occupation;
            particles.occupation[static_cast<size_t>(i)] = std::max(0.0, new_occupation);
        }
        residual_raw_ev += (energy_after - energy_before) * cell_energy_weight;
    }

    // Grid cells are constrained to one material, so accumulate the physical
    // residual directly rather than applying one global material weight.
    return residual_raw_ev;
}

}  // namespace phonomc
