#pragma once

#include "SimulationDomain.h"

#include <array>
#include <limits>
#include <stdexcept>

namespace phonomc {
// Shared by initialization, transport and injection. The fallback preserves
// first-center tie breaking and the original squared-distance evaluation order.
inline int nearest_grid_index(const SimulationDomain& geometry,
                              const std::array<double, 3>& position) {
    if (geometry.fast_grid_index_enabled()) {
        const int idx = geometry.fast_grid_index(position);
        if (idx >= 0) {
            return idx;
        }
    }
    const auto& centers = geometry.grid_centers();
    if (centers.empty()) {
        throw std::logic_error("Grid lookup requires a nonempty domain grid.");
    }
    int best = 0;
    double best_d2 = std::numeric_limits<double>::max();
    for (int i = 0; i < static_cast<int>(centers.size()); ++i) {
        const std::array<double, 3> delta{
            position[0] - centers[i][0], position[1] - centers[i][1], position[2] - centers[i][2]};
        const double d2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];
        if (d2 < best_d2) {
            best_d2 = d2;
            best = i;
        }
    }
    return best;
}

} // namespace phonomc
