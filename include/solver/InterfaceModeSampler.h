#pragma once
#include <array>
#include <random>
#include <vector>

class PhononMaterial;
namespace phonomc {
// Immutable directional flux tables after configure(). Sampling changes only
// the caller's RNG, so separate thread-local RNGs can share these tables.
class InterfaceModeSampler {
public:
    void configure(const std::vector<const PhononMaterial*>& materials, int axis);
    int sample(int material_id, int direction, double omega, double tolerance,
               std::mt19937_64& rng, double& flux_sum) const;
    // Shared half-open frequency bands, independent of incident material/mode.
    int sample_band(int material_id, int direction, double omega, double band_width,
                    std::mt19937_64& rng, double& flux_sum) const;
private:
    struct Bank {
        std::vector<int> active_indices;
        std::vector<double> angular_frequencies;
        std::vector<double> flux_prefix;
    };
    std::vector<std::array<Bank, 2>> banks_;
    int draw(const Bank& bank, size_t begin, size_t end,
             std::mt19937_64& rng, double& flux_sum) const;
};
}  // namespace phonomc
