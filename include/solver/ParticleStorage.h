#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <tuple>
#include <vector>

namespace phonomc {

using Vec3 = std::array<double, 3>;

// AoS only at insertion boundaries; transport retains contiguous SoA arrays.
struct ParticleRecord {
    std::array<int, 2> mode = {0, 0};
    int material_id = 0;
    Vec3 position = {};
    Vec3 velocity = {};
    Vec3 collision_position = {};
    double collision_time = std::numeric_limits<double>::infinity();
    double temperature = 300.0;
    double omega = 0.0;
    double occupation = 0.0;
    double energy = 0.0;
    int grid_id = 0;
    int collision_facet = -1;
    char collision_condition = 'R';
    std::uint8_t alive = 1;
    std::uint8_t collision_failed = 0;
};

// Owns all per-particle fields, including the collision cache. Call reset,
// append or compact_alive for structural changes, never resize an individual
// field. Kernels may update existing elements directly without accessor overhead.
// Structural operations must run outside parallel particle loops.
class ParticleStorage {
public:
    int size() const { return static_cast<int>(modes.size()); }
    bool aligned() const;
    void reset(int count);
    int append(const ParticleRecord& particle);
    int compact_alive();

    // Fixed quadrature volumes by material. Empty retains the single-material
    // API default. Multimaterial sampling uses equal physical state weights.
    std::vector<double> material_spatial_volumes_a3;
    double spatial_volume(int material, double fallback) const {
        return material_spatial_volumes_a3.empty() ? fallback :
            material_spatial_volumes_a3.at(static_cast<size_t>(material));
    }

    std::vector<std::array<int, 2>> modes;
    std::vector<int> material_ids;
    std::vector<Vec3> positions;
    std::vector<Vec3> velocities;
    std::vector<Vec3> collision_positions;
    std::vector<double> collision_times;
    std::vector<double> temperatures;
    std::vector<double> omega;
    std::vector<double> occupation;
    std::vector<double> energies;
    std::vector<int> grid_ids;
    std::vector<int> collision_facets;
    std::vector<char> collision_conditions;
    std::vector<std::uint8_t> alive;
    std::vector<std::uint8_t> collision_failed;

private:
    auto payload_arrays() {
        return std::tie(modes, material_ids, positions, velocities,
                        collision_positions, collision_times, temperatures,
                        omega, occupation, energies, grid_ids, collision_facets,
                        collision_conditions, collision_failed);
    }

    auto arrays() {
        return std::tuple_cat(payload_arrays(), std::tie(alive));
    }
};

}  // namespace phonomc
