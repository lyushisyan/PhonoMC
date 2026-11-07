#pragma once
#include <array>
#include <limits>
#include <vector>

class SimulationDomain;

namespace phonomc {
struct MaterialCrossing {
    double time_ps = std::numeric_limits<double>::infinity();
    int destination = -1;
};

// Immutable spatial layout. No particles, material physics, RNG or solver state.
class LayeredMaterialLayout {
public:
    LayeredMaterialLayout() = default;
    LayeredMaterialLayout(int material_count, int axis, std::vector<double> planes);
    int material_at(const std::array<double, 3>& position) const;
    MaterialCrossing next_crossing(const std::array<double, 3>& position,
                                   const std::array<double, 3>& velocity,
                                   int current_material) const;
    // Also checks alignment, occupied layers and compatible boundary faces.
    std::vector<int> validated_cell_materials(const SimulationDomain& domain,
                                              std::array<int, 3> grid_shape) const;

    int axis() const { return axis_; }

    bool has_interfaces() const { return material_count_ > 1 && !planes_.empty(); }

private:
    int material_count_ = 1;
    int axis_ = 0;
    std::vector<double> planes_;
};
} // namespace phonomc
