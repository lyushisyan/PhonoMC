#include "solver/ReservoirSampler.h"
#include "PhononMaterial.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace phonomc {
namespace {
double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
}
void ReservoirSampler::configure(const std::vector<ReservoirFacetSpec>& facets,
    const std::vector<const PhononMaterial*>& materials, int mesh_facet_count,
    double particle_density, double time_step, const std::vector<double>& material_particle_densities) {
    if (mesh_facet_count < 0 || !std::isfinite(particle_density) || particle_density < 0 ||
        !std::isfinite(time_step) || time_step < 0)
        throw std::invalid_argument("Invalid reservoir discretization.");
    if (!material_particle_densities.empty() && material_particle_densities.size()!=materials.size())
        throw std::invalid_argument("Reservoir material density count mismatch.");
    for (double density:material_particle_densities)
        if (!std::isfinite(density) || density<=0)
            throw std::invalid_argument("Invalid reservoir material density.");
    ReservoirSampler next;
    next.time_step_ = time_step;
    next.facet_to_reservoir_.assign(static_cast<size_t>(mesh_facet_count), -1);
    for (auto spec : facets) {
        if (spec.facet < 0 || spec.facet >= mesh_facet_count || spec.material < 0 ||
            spec.material >= static_cast<int>(materials.size()) || materials[spec.material] == nullptr)
            throw std::invalid_argument("Invalid reservoir facet/material index.");
        if (next.facet_to_reservoir_[spec.facet] >= 0)
            throw std::invalid_argument("Duplicate reservoir facet.");
        if (!std::isfinite(spec.temperature) || spec.temperature < 0 ||
            !std::isfinite(spec.area) || spec.area < 0)
            throw std::invalid_argument("Invalid reservoir temperature/area.");
        for (double x : spec.outward_normal)
            if (!std::isfinite(x)) throw std::invalid_argument("Non-finite reservoir normal.");
        if (std::abs(std::sqrt(dot(spec.outward_normal, spec.outward_normal)) - 1.0) > 1e-8)
            throw std::invalid_argument("Reservoir normal must have unit length.");
        const auto& phonon = *materials[spec.material];
        const auto& modes = phonon.active_mode_list();
        if (modes.empty()) throw std::invalid_argument("Reservoir material has no active modes.");
        spec.area = std::max(1e-12, spec.area);
        Entry entry;
        entry.facet = spec;
        entry.cdf.resize(modes.size());
        const double n_active = std::max(1.0, static_cast<double>(phonon.active_mode_count()));
        const double density = material_particle_densities.empty() ? particle_density :
            material_particle_densities[spec.material];
        const double thickness = n_active / std::max(1e-18, density * spec.area);
        double expected = 0;
        for (size_t m = 0; m < modes.size(); ++m) {
            const double parallel = -dot(spec.outward_normal, phonon.mode_group_velocity(modes[m]));
            expected += std::max(0.0, parallel * time_step / std::max(1e-18, thickness));
            entry.cdf[m] = expected;
        }
        if (!std::isfinite(expected) || expected > static_cast<double>(std::numeric_limits<int>::max()))
            throw std::overflow_error("Reservoir initial emission count exceeds int capacity.");
        if (expected > 0) for (double& value : entry.cdf) value /= expected;
        entry.initial_count = std::max(0, static_cast<int>(std::llround(expected)));
        entry.expected_count = expected;
        next.facet_to_reservoir_[spec.facet] = next.size();
        next.entries_.push_back(std::move(entry));
    }
    *this = std::move(next);
}
int ReservoirSampler::reservoir_index(int facet) const {
    return facet >= 0 && facet < static_cast<int>(facet_to_reservoir_.size())
        ? facet_to_reservoir_[facet] : -1;
}
const ReservoirFacetSpec& ReservoirSampler::facet(int reservoir) const {
    return entries_.at(static_cast<size_t>(reservoir)).facet;
}
int ReservoirSampler::initial_emission_count(int reservoir) const {
    return entries_.at(static_cast<size_t>(reservoir)).initial_count;
}
int ReservoirSampler::sample_flux_emission_count(int reservoir, std::mt19937_64& rng) const {
    const double expected = entries_.at(static_cast<size_t>(reservoir)).expected_count;
    const int count = static_cast<int>(std::floor(expected));
    const double fraction = expected - count;
    return count + (fraction > 0.0 && std::uniform_real_distribution<double>(0.0, 1.0)(rng) < fraction ? 1 : 0);
}
std::vector<ReservoirEmission> ReservoirSampler::sample(int reservoir, int count, std::mt19937_64& rng) const {
    const auto& roulette = entries_.at(static_cast<size_t>(reservoir)).cdf;
    std::vector<ReservoirEmission> emissions;
    if (count <= 0) return emissions;
    emissions.reserve(static_cast<size_t>(count));
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    if (roulette.back() > 0) {
        for (int c = 0; c < count; ++c) {
            const double rr = U01(rng);
            const auto it = std::lower_bound(roulette.begin(), roulette.end(), rr);
            const size_t m = it == roulette.end() ? roulette.size()-1 : static_cast<size_t>(it-roulette.begin());
            emissions.push_back({m, time_step_ * U01(rng)});
        }
    } else {
        // Preserve the legacy zero-flux fallback, including its random draw order.
        std::uniform_int_distribution<int> Uidx(0, static_cast<int>(roulette.size())-1);
        for (int c = 0; c < count; ++c) {
            const auto m = static_cast<size_t>(Uidx(rng));
            emissions.push_back({m, time_step_ * U01(rng)});
        }
    }
    return emissions;
}
} // namespace phonomc
