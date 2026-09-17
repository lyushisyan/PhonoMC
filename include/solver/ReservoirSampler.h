#pragma once
#include <array>
#include <random>
#include <vector>

class PhononMaterial;
namespace phonomc {
struct ReservoirFacetSpec {
    int facet = -1, material = 0;
    double temperature = 300, area = 0;
    std::array<double, 3> outward_normal {};
};
struct ReservoirEmission {
    size_t active_index = 0;
    double time_ps = 0;
};

// Precomputed incoming-flux CDFs. No geometry, particle storage or owned RNG.
class ReservoirSampler {
public:
    void configure(const std::vector<ReservoirFacetSpec>& facets,
                   const std::vector<const PhononMaterial*>& materials,
                   int mesh_facet_count, double particle_density, double time_step,
                   const std::vector<double>& material_particle_densities = {});
    int size() const { return static_cast<int>(entries_.size()); }
    int reservoir_index(int facet) const;
    const ReservoirFacetSpec& facet(int reservoir) const;
    int initial_emission_count(int reservoir) const;
    // Unbiased fixed-flux count, independent of the previous absorbed population.
    int sample_flux_emission_count(int reservoir, std::mt19937_64& rng) const;
    // Draw all mode/time pairs before the caller samples any surface positions.
    std::vector<ReservoirEmission> sample(int reservoir, int count, std::mt19937_64& rng) const;
private:
    struct Entry {
        ReservoirFacetSpec facet;
        std::vector<double> cdf;
        int initial_count = 0;
        double expected_count = 0;
    };
    std::vector<Entry> entries_;
    std::vector<int> facet_to_reservoir_;
    double time_step_ = 0;
};
} // namespace phonomc
