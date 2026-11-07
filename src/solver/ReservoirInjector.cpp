#include "solver/ReservoirInjector.h"

#include "PhononMaterial.h"
#include "SimulationDomain.h"
#include "solver/BoundaryLedger.h"
#include "solver/GridLocation.h"
#include "solver/ParticleStorage.h"
#include "solver/ReservoirSampler.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace phonomc {
namespace {
Vec3 add(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 mul(const Vec3& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}
} // namespace

ReservoirInjector::ReservoirInjector(ReservoirInjectionOptions options,
                                     const SimulationDomain& geometry,
                                     const std::vector<const PhononMaterial*>& materials,
                                     std::vector<double> energy_weights,
                                     const ReservoirSampler& sampler)
    : options_(options), geometry_(geometry), materials_(materials),
      energy_weights_(std::move(energy_weights)), sampler_(sampler) {
    if (!std::isfinite(options_.reference_temperature) || options_.reference_temperature < 0 ||
        !std::isfinite(options_.push_distance_a) || options_.push_distance_a < 0 ||
        materials_.empty() || energy_weights_.size() != materials_.size() ||
        (options_.generation != ReservoirGeneration::OneToOne &&
         options_.generation != ReservoirGeneration::FixedIncomingFlux)) {
        throw std::invalid_argument("Invalid reservoir injection options or material weights.");
    }
    for (size_t i = 0; i < materials_.size(); ++i) {
        if (!materials_[i] || !std::isfinite(energy_weights_[i]) || energy_weights_[i] <= 0) {
            throw std::invalid_argument("Invalid reservoir injection material or energy weight.");
        }
    }
    for (int r = 0; r < sampler_.size(); ++r) {
        const auto& entry = sampler_.facet(r);
        if (entry.material < 0 || entry.material >= static_cast<int>(materials_.size()) ||
            entry.facet < 0 || entry.facet >= geometry_.mesh().facet_count()) {
            throw std::invalid_argument(
                "Reservoir sampler does not match injection materials/geometry.");
        }
    }
}

std::vector<std::pair<int, double>> ReservoirInjector::inject(ParticleStorage& particles,
                                                              BoundaryLedger& ledger,
                                                              std::mt19937_64& rng) const {
    std::vector<std::pair<int, double>> inserted;
    if (sampler_.size() <= 0) {
        return inserted;
    }

    const auto& mesh = geometry_.mesh();

    for (int reservoir = 0; reservoir < sampler_.size(); ++reservoir) {
        const auto& entry = sampler_.facet(reservoir);
        const int material_id = entry.material;
        const PhononMaterial& phonon = *materials_[material_id];
        const auto& modes = phonon.active_mode_list();
        // Multimaterial DMM needs prescribed incoming phase-space density,
        // independent of the material-dependent outgoing carrier population.
        const int emission_count = (options_.generation == ReservoirGeneration::FixedIncomingFlux)
                                       ? sampler_.sample_flux_emission_count(reservoir, rng)
                                       : ledger.emission_count(reservoir);
        const auto emissions = sampler_.sample(reservoir, emission_count, rng);
        if (emissions.empty()) {
            continue;
        }
        // Keep position sampling after the complete mode/time batch.
        std::vector<Vec3> surface_positions = mesh.sample_surface_points(
            static_cast<int>(emissions.size()), std::vector<int>{entry.facet}, rng);
        const Vec3 outward_normal = entry.outward_normal;
        const double reservoir_temperature = entry.temperature;

        for (size_t k = 0; k < emissions.size(); ++k) {
            const std::array<int, 2> mode = modes[emissions[k].active_index];
            Vec3 velocity = phonon.mode_group_velocity(mode);
            if (dot(velocity, outward_normal) > 0.0) {
                if (options_.conservative) {
                    throw std::runtime_error("Full-matrix reservoir sampled an outward mode.");
                }
                velocity = mul(velocity, -1.0);
            }
            const double speed = std::max(1e-12, norm(velocity));
            Vec3 position =
                add(surface_positions[k], mul(velocity, options_.push_distance_a / speed));

            const double occupation =
                options_.linearized
                    ? phonon.linearized_bose_occupation(
                          reservoir_temperature, options_.reference_temperature, mode)
                    : phonon.bose_occupation(reservoir_temperature, mode);
            phonomc::ParticleRecord particle;
            particle.mode = mode;
            particle.material_id = material_id;
            particle.position = position;
            particle.velocity = velocity;
            particle.collision_position = position;
            particle.temperature = reservoir_temperature;
            particle.omega = phonon.mode_angular_frequency(mode);
            particle.occupation = occupation;
            particle.grid_id = phonomc::nearest_grid_index(geometry_, position);
            const int particle_index = particles.append(particle);
            inserted.push_back({particle_index, emissions[k].time_ps});

            constexpr double kHbarEvPs = 6.582119569e-4;
            const double equilibrium_occupation =
                phonon.bose_occupation(options_.reference_temperature, mode);
            const double raw_deviation_ev = kHbarEvPs * phonon.mode_angular_frequency(mode) *
                                            (occupation - equilibrium_occupation);
            ledger.record_injection(raw_deviation_ev * energy_weights_[material_id]);
        }
    }
    return inserted;
}
} // namespace phonomc
