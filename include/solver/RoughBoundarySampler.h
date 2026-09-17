#pragma once
#include <array>
#include <ostream>
#include <random>
#include <vector>

class PhononMaterial;
namespace phonomc {
struct RoughFacetSpec {
    bool enabled = false;
    std::array<double, 3> outward_normal {};
    double roughness = 0;
};
struct RoughSampleEvents {
    long long events = 0, specular = 0, diffuse = 0;
    long long residual_window = 0, residual_fallback = 0;
    long long missing_data = 0, missing_spec_match = 0, outgoing_pool = 0, global_random = 0;
};
struct RoughReflection {
    std::array<int, 2> mode {};
    double occupation = 0;
    RoughSampleEvents diagnostics;
};

// Owns immutable lookup tables, not geometry, particles, RNGs or atomic counters.
// The material passed to sample() must match the original table's mode ordering.
class RoughBoundarySampler {
public:
    void configure(const std::vector<RoughFacetSpec>& facets,
                   const std::vector<const PhononMaterial*>& materials,
                   std::ostream* progress = nullptr, bool require_discrete_reflection = false,
                   double equilibrium_temperature = 0.0);
    int table_index(int material, int facet) const;
    RoughReflection sample(const PhononMaterial& material, int table,
                           const std::array<int, 2>& mode, double occupation,
                           double reference_temperature, std::mt19937_64& rng,
                           bool signed_weights = false) const;
private:
    using Vec3 = std::array<double, 3>;
    struct RoughFacetData {
        int facet = -1, material_id = 0;
        bool strict_elastic = false;
        std::vector<double> specularity;
        std::vector<int> spec_match_active, outgoing_active, outgoing_sorted_active;
        std::vector<double> outgoing_sorted_omega, outgoing_sorted_residual_flux_prefix;
        std::vector<int> diffuse_begin, diffuse_end, diffuse_roulette_active;
        std::vector<double> diffuse_roulette_cdf;
    };
    void build(const std::vector<RoughFacetSpec>& facets,
               const std::vector<const PhononMaterial*>& materials, std::ostream* progress,
               bool require_discrete_reflection, double equilibrium_temperature);
    int sample_diffuse(int table, int incoming, int* source,
                       std::mt19937_64& rng, RoughSampleEvents& events) const;
    int facet_count_ = 0;
    std::vector<int> active_counts_, facet_to_rough_data_;
    std::vector<RoughFacetData> rough_boundary_data_;
};
}  // namespace phonomc
