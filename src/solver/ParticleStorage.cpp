#include "solver/ParticleStorage.h"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <utility>

namespace phonomc {

bool ParticleStorage::aligned() const {
    return material_ids.size() == modes.size() &&
           positions.size() == modes.size() &&
           velocities.size() == modes.size() &&
           collision_positions.size() == modes.size() &&
           collision_times.size() == modes.size() &&
           temperatures.size() == modes.size() &&
           omega.size() == modes.size() &&
           occupation.size() == modes.size() &&
           energies.size() == modes.size() &&
           grid_ids.size() == modes.size() &&
           collision_facets.size() == modes.size() &&
           collision_conditions.size() == modes.size() &&
           alive.size() == modes.size() &&
           collision_failed.size() == modes.size();
}

void ParticleStorage::reset(int count) {
    if (count < 0) throw std::invalid_argument("Particle count must be nonnegative.");
    // Build separately so allocation failure leaves the existing storage valid.
    ParticleStorage fresh;
    const auto n = static_cast<std::size_t>(count);
    std::apply([&](auto&... values) { (values.resize(n), ...); }, fresh.arrays());
    std::fill(fresh.collision_times.begin(), fresh.collision_times.end(),
              std::numeric_limits<double>::infinity());
    std::fill(fresh.temperatures.begin(), fresh.temperatures.end(), 300.0);
    std::fill(fresh.collision_facets.begin(), fresh.collision_facets.end(), -1);
    std::fill(fresh.collision_conditions.begin(), fresh.collision_conditions.end(), 'R');
    std::fill(fresh.alive.begin(), fresh.alive.end(), std::uint8_t{1});
    *this = std::move(fresh);
}

int ParticleStorage::append(const ParticleRecord& p) {
    assert(aligned());
    const auto n = modes.size();
    const auto limit = static_cast<std::size_t>(std::numeric_limits<int>::max());
    if (n >= limit) throw std::length_error("Particle index exceeds int range.");
    // Reserve every field before changing any size. The scalar/array pushes
    // below then cannot allocate or throw; all fields stay aligned on failure.
    const auto capacity = std::max(n + 1, std::min(limit, n + n / 2));
    std::apply([&](auto&... values) {
        ((values.capacity() < n + 1 ? values.reserve(capacity) : void()), ...);
    }, arrays());
    modes.push_back(p.mode);
    material_ids.push_back(p.material_id);
    positions.push_back(p.position);
    velocities.push_back(p.velocity);
    collision_positions.push_back(p.collision_position);
    collision_times.push_back(p.collision_time);
    temperatures.push_back(p.temperature);
    omega.push_back(p.omega);
    occupation.push_back(p.occupation);
    energies.push_back(p.energy);
    grid_ids.push_back(p.grid_id);
    collision_facets.push_back(p.collision_facet);
    collision_conditions.push_back(p.collision_condition);
    alive.push_back(p.alive);
    collision_failed.push_back(p.collision_failed);
    return static_cast<int>(n);
}

int ParticleStorage::compact_alive() {
    assert(aligned());
    const auto kept = static_cast<std::size_t>(
        std::count_if(alive.begin(), alive.end(), [](std::uint8_t a) { return a != 0; }));
    const int removed = static_cast<int>(alive.size() - kept);
    if (removed == 0) return 0;
    const auto compact = [&](auto& values) {
        // Keep the mask unchanged until every payload has been compacted.
        std::size_t write = 0;
        for (std::size_t read = 0; read < alive.size(); ++read) {
            if (alive[read] != 0) {
                if (write != read) values[write] = std::move(values[read]);
                ++write;
            }
        }
        values.resize(kept);
    };
    std::apply([&](auto&... values) { (compact(values), ...); }, payload_arrays());
    alive.assign(kept, std::uint8_t{1});
    return removed;
}

}  // namespace phonomc
