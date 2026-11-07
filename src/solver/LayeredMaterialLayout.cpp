#include "solver/LayeredMaterialLayout.h"

#include "SimulationDomain.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace phonomc {
LayeredMaterialLayout::LayeredMaterialLayout(int count, int axis, std::vector<double> planes)
    : material_count_(count), axis_(axis), planes_(std::move(planes)) {
    if (count < 1 || axis < 0 || axis > 2 || planes_.size() != static_cast<size_t>(count - 1)) {
        throw std::invalid_argument("Invalid layered material count, axis or interface count.");
    }
    for (size_t i = 0; i < planes_.size(); ++i) {
        if (!std::isfinite(planes_[i]) || (i && planes_[i] <= planes_[i - 1])) {
            throw std::invalid_argument(
                "Material interface positions must be finite and strictly increasing.");
        }
    }
}

int LayeredMaterialLayout::material_at(const std::array<double, 3>& position) const {
    if (!has_interfaces()) {
        return 0;
    }
    const double coordinate = position[static_cast<size_t>(axis_)];
    return static_cast<int>(std::upper_bound(planes_.begin(), planes_.end(), coordinate) -
                            planes_.begin());
}

MaterialCrossing LayeredMaterialLayout::next_crossing(const std::array<double, 3>& position,
                                                      const std::array<double, 3>& velocities,
                                                      int current) const {
    MaterialCrossing crossing;
    if (material_count_ <= 1) {
        return crossing;
    }
    const double velocity = velocities[static_cast<size_t>(axis_)];
    if (std::abs(velocity) <= 1e-18) {
        return crossing;
    }
    double plane = 0.0;
    if (velocity > 0.0 && current + 1 < material_count_) {
        plane = planes_[static_cast<size_t>(current)];
        crossing.destination = current + 1;
    } else if (velocity < 0.0 && current > 0) {
        plane = planes_[static_cast<size_t>(current - 1)];
        crossing.destination = current - 1;
    } else {
        return crossing;
    }
    const double time = (plane - position[static_cast<size_t>(axis_)]) / velocity;
    if (!std::isfinite(time) || time <= 1e-14) {
        return {};
    }
    crossing.time_ps = time;
    return crossing;
}

std::vector<int>
LayeredMaterialLayout::validated_cell_materials(const SimulationDomain& geometry,
                                                std::array<int, 3> grid_shape) const {
    std::vector<int> cell_materials(geometry.grid_centers().size(), 0);
    for (size_t grid = 0; grid < geometry.grid_centers().size(); ++grid) {
        cell_materials[grid] = material_at(geometry.grid_centers()[grid]);
    }
    if (static_cast<size_t>(material_count_) <= 1) {
        return cell_materials;
    }
    const int cells = (axis_ == 0) ? grid_shape[0] : (axis_ == 1) ? grid_shape[1] : grid_shape[2];
    const double lo = geometry.bounds_min()[static_cast<size_t>(axis_)];
    const double hi = geometry.bounds_max()[static_cast<size_t>(axis_)];
    const double cell_width = (hi - lo) / static_cast<double>(cells);
    const double tolerance = std::max(1e-10, 1e-10 * std::abs(hi - lo));
    for (double plane : planes_) {
        if (!std::isfinite(plane) || !(plane > lo) || !(plane < hi)) {
            throw std::runtime_error(
                "Each material interface must lie strictly inside the loaded geometry bounds.");
        }
        const double cell_coordinate = (plane - lo) / cell_width;
        if (std::abs(cell_coordinate - std::round(cell_coordinate)) > tolerance / cell_width) {
            throw std::runtime_error("Each materials.interface_positions value must coincide with "
                                     "a grid-cell boundary.");
        }
    }
    std::vector<bool> occupied(static_cast<size_t>(material_count_), false);
    for (int id : cell_materials) {
        occupied.at(static_cast<size_t>(id)) = true;
    }
    if (std::find(occupied.begin(), occupied.end(), false) != occupied.end()) {
        throw std::runtime_error(
            "Every material layer must contain at least one active geometry cell.");
    }
    const auto& normals = geometry.mesh().facet_normals();
    const int axis = axis_;
    for (int facet = 0; facet < geometry.mesh().facet_count(); ++facet) {
        const char condition = geometry.facet_boundary_condition(facet);
        if (condition == 'P' &&
            std::abs(normals[static_cast<size_t>(facet)][static_cast<size_t>(axis)]) > 0.99) {
            throw std::runtime_error(
                "Periodic boundaries normal to materials.interface_axis are not supported, "
                "because they would join different end materials.");
        }
        if (condition == 'T' &&
            std::abs(normals[static_cast<size_t>(facet)][static_cast<size_t>(axis)]) < 0.99) {
            throw std::runtime_error(
                "In the layered-material model, thermal-reservoir facets must be normal to "
                "materials.interface_axis so each reservoir touches one material.");
        }
    }
    return cell_materials;
}
} // namespace phonomc
