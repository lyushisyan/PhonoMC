#pragma once

#include "PhononMaterial.h"
#include <stdexcept>
#include <vector>

namespace phonomc {
// Immutable after construction: tied to one material and one fixed background.
// Dense q/branch indexing also covers inactive modes without an active-map lookup.
class BackgroundCache {
public:
    BackgroundCache(const PhononMaterial& material, double temperature, double particle_volume,
                    bool exact_reference = false)
        : branches_(material.branch_count()),
          energy_weight_(static_cast<double>(material.active_mode_count()) * particle_volume /
                         material.energy_density_normalization()),
          energy_density_(exact_reference ? material.crystal_energy_density(temperature)
                                          : material.energy_density_from_temperature(temperature)), temperature_(temperature) {
        occupation_.reserve(static_cast<std::size_t>(material.qpoint_count()) * branches_);
        for (int q = 0; q < material.qpoint_count(); ++q)
            for (int b = 0; b < branches_; ++b) {
                occupation_.push_back(material.bose_occupation(temperature, {q, b}));
                heat_capacity_density_ += material.mode_heat_capacity(temperature, {q, b}) /
                    material.energy_density_normalization();
            }
    }
    double occupation(const PhononMaterial::Mode& mode) const {
        return occupation_.at(static_cast<std::size_t>(mode[0]) * branches_ + mode[1]);
    }
    bool contains_mode(const PhononMaterial::Mode& mode) const noexcept {
        return branches_ > 0 && mode[0] >= 0 && mode[1] >= 0 && mode[1] < branches_ &&
            static_cast<std::size_t>(mode[0]) < occupation_.size() / static_cast<std::size_t>(branches_);
    }
    double energy_weight() const { return energy_weight_; }
    double energy_density() const { return energy_density_; }
    double temperature() const { return temperature_; }
    double heat_capacity_density() const { return heat_capacity_density_; }
private:
    int branches_;
    double energy_weight_, energy_density_;
    double temperature_, heat_capacity_density_ = 0.0;
    std::vector<double> occupation_;
};
}  // namespace phonomc
