#include "solver/PrescribedHeatSource.h"

#include "SimulationConfig.h"
#include "SimulationDomain.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <sstream>

namespace phonomc {
namespace {
using Vec3 = std::array<double, 3>;
constexpr double kEvPerWattPs = 6.241509074e6;
}

void PrescribedHeatSource::configure(const SimulationConfig& args_, const SimulationDomain& geometry) {
    *this = PrescribedHeatSource{};
    if (!args_.heat_source_enabled) return;
    validate_heat_source_config(args_);
    const auto profile = args_.heat_source_profile;
    const double density = args_.heat_source_total_power.has_value() ? 1.0 : args_.heat_source_power_density;
    std::vector<double> weights;

    // Heat-source min/max/center/sigma are absolute lengths in the input file
    // (nm) and have already been converted to internal Angstrom by the
    // configuration loader.
    const auto& centers = geometry.grid_centers();
    positions_=centers;
    weights.assign(centers.size(), 0.0);
    int selected = 0;

    if (profile == HeatSourceProfile::Uniform) {
        if (args_.heat_source_min.size() != 3 || args_.heat_source_max.size() != 3) {
            throw std::runtime_error("uniform heat source requires min/max 3D vectors.");
        }
        Vec3 hs_min = {
            args_.heat_source_min[0], args_.heat_source_min[1], args_.heat_source_min[2]
        };
        Vec3 hs_max = {
            args_.heat_source_max[0], args_.heat_source_max[1], args_.heat_source_max[2]
        };
        for (int k = 0; k < 3; ++k) {
            if (hs_min[k] > hs_max[k]) {
                std::swap(hs_min[k], hs_max[k]);
            }
        }
        for (size_t i = 0; i < centers.size(); ++i) {
            const Vec3 c = centers[i];
            const bool inside =
                c[0] >= hs_min[0] && c[0] <= hs_max[0] &&
                c[1] >= hs_min[1] && c[1] <= hs_max[1] &&
                c[2] >= hs_min[2] && c[2] <= hs_max[2];
            if (inside) {
                weights[i] = 1.0;
                ++selected;
            }
        }
        if (selected == 0) {
            throw std::runtime_error("uniform heat source region does not include any grid center.");
        }
    } else {
        if ((args_.heat_source_center.size() != 3 && args_.heat_source_centers.empty()) ||
            args_.heat_source_sigma.size() != 3 ||
            (!args_.heat_source_centers.empty() && (args_.heat_source_centers.size() % 3) != 0)) {
            throw std::runtime_error("gaussian heat source requires center=[x,y,z] or centers=[[x,y,z],...] and sigma=[sx,sy,sz].");
        }
        std::vector<Vec3> hs_centers;
        if (!args_.heat_source_centers.empty()) {
            hs_centers.reserve(args_.heat_source_centers.size() / 3);
            for (size_t j = 0; j + 2 < args_.heat_source_centers.size(); j += 3) {
                hs_centers.push_back({
                    args_.heat_source_centers[j],
                    args_.heat_source_centers[j + 1],
                    args_.heat_source_centers[j + 2]
                });
            }
        } else {
            hs_centers.push_back({
                args_.heat_source_center[0],
                args_.heat_source_center[1],
                args_.heat_source_center[2]
            });
        }
        std::vector<double> source_peak_scale(hs_centers.size(), 1.0);
        if (!args_.heat_source_power_densities.empty()) {
            if (args_.heat_source_power_densities.size() != hs_centers.size()) {
                throw std::runtime_error("heat_source.power_densities size does not match heat_source.centers.");
            }
            for (size_t j = 0; j < hs_centers.size(); ++j) {
                const double qj = args_.heat_source_power_densities[j];
                if (!(std::isfinite(qj) && qj > 0.0)) {
                    throw std::runtime_error("heat_source.power_densities must be finite positive values.");
                }
                source_peak_scale[j] = qj / density;
            }
        }
        Vec3 hs_sigma {1.0, 1.0, 1.0};
        Vec3 hs_min {0.0, 0.0, 0.0};
        Vec3 hs_max {0.0, 0.0, 0.0};
        const bool has_clip_box = args_.heat_source_min.size() == 3 && args_.heat_source_max.size() == 3;
        for (int k = 0; k < 3; ++k) {
            if (has_clip_box) {
                hs_min[k] = args_.heat_source_min[k];
                hs_max[k] = args_.heat_source_max[k];
                if (hs_min[k] > hs_max[k]) {
                    std::swap(hs_min[k], hs_max[k]);
                }
            }
        }
        std::array<bool, 3> gaussian_axis_enabled {true, true, true};
        std::array<double, 3> uniform_half_width {
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()
        };
        if (!args_.heat_source_half_width.empty()) {
            if (args_.heat_source_half_width.size() != 3) {
                throw std::runtime_error("heat_source.half_width must contain 3 values.");
            }
            for (int k = 0; k < 3; ++k) {
                if (std::isfinite(args_.heat_source_half_width[k]) && args_.heat_source_half_width[k] > 0.0) {
                    uniform_half_width[k] = args_.heat_source_half_width[k];
                }
            }
        }
        for (int k = 0; k < 3; ++k) {
            if (std::isfinite(args_.heat_source_sigma[k]) && args_.heat_source_sigma[k] > 0.0) {
                hs_sigma[k] = std::max(1e-12, args_.heat_source_sigma[k]);
                gaussian_axis_enabled[k] = true;
            } else {
                // sigma<=0 means this axis is uniform (no Gaussian variation).
                hs_sigma[k] = 1.0;
                gaussian_axis_enabled[k] = false;
            }
        }
        for (size_t i = 0; i < centers.size(); ++i) {
            const Vec3 c = centers[i];
            bool outside_cutoff = false;
            if (has_clip_box) {
                for (int k = 0; k < 3; ++k) {
                    if (c[k] < hs_min[k] || c[k] > hs_max[k]) {
                        outside_cutoff = true;
                        break;
                    }
                }
            }
            if (outside_cutoff) {
                weights[i] = 0.0;
                continue;
            }
            double wsum = 0.0;
            for (size_t src = 0; src < hs_centers.size(); ++src) {
                double r2 = 0.0;
                bool outside_source = false;
                for (int k = 0; k < 3; ++k) {
                    if (gaussian_axis_enabled[k]) {
                        const double u = std::abs(c[k] - hs_centers[src][k]) / hs_sigma[k];
                        if (u > 2.0) {
                            outside_source = true;
                            break;
                        }
                        r2 += u * u;
                    } else if (std::isfinite(uniform_half_width[k]) &&
                               std::abs(c[k] - hs_centers[src][k]) > uniform_half_width[k]) {
                        outside_source = true;
                        break;
                    }
                }
                if (!outside_source) {
                    wsum += source_peak_scale[src] * std::exp(-0.5 * r2);
                }
            }
            if (wsum > 0.0) {
                weights[i] = wsum;
                ++selected;
            }
        }
        if (selected == 0) {
            throw std::runtime_error("gaussian heat source cannot be applied because grid is empty.");
        }
    }

    const auto& volumes = geometry.grid_volumes();
    double weighted_volume = 0.0;
    for (size_t cell = 0; cell < weights.size(); ++cell)
        weighted_volume += weights[cell] * volumes[cell];
    if (!(weighted_volume > 0.0) || !std::isfinite(weighted_volume))
        throw std::runtime_error("Heat source has no finite positive sampled volume; refine the grid or source region.");
    rates_.resize(weights.size());
    for (size_t cell = 0; cell < weights.size(); ++cell) {
        rates_[cell] = args_.heat_source_total_power.has_value()
            ? *args_.heat_source_total_power * kEvPerWattPs * (weights[cell] * volumes[cell] / weighted_volume)
            : density * 6.241509074e-24 * weights[cell] * volumes[cell];
        if (!std::isfinite(rates_[cell]) || rates_[cell] < 0.0)
            throw std::runtime_error("Non-finite heat-source cell power.");
    }
    const double total_rate = std::accumulate(rates_.begin(), rates_.end(), 0.0);
    if (!std::isfinite(total_rate) || !(total_rate > 0.0))
        throw std::runtime_error("Total heat-source rate must be finite and representably positive.");
    prescribed_.assign(rates_.size(), 0.0);
    deposited_.assign(rates_.size(), 0.0);
    pending_.assign(rates_.size(), 0.0);
    enabled_ = true;
    // Solver initialization previously selected fixed two-decimal formatting.
    // Format power independently so micro/nanowatt inputs are not printed as zero.
    std::ostringstream power_text;
    power_text << std::scientific << std::setprecision(8) << base_power_w();
    std::cout << "Prescribed lattice heat source: cells=" << selected
              << ", profile=" << to_string(profile)
              << ", base_power=" << power_text.str() << " W"
              << ", power_input=" << (args_.heat_source_total_power ? "total_power" : "power_density")
              << ", empty_cell_policy=" << (args_.collision_model == CollisionModel::FullMatrix ? "create_modal_carriers" : "defer_locally")
              << ", electron_transport=false\n";
}

void PrescribedHeatSource::accrue(double integrated_time_ps) {
    if (!std::isfinite(integrated_time_ps) || integrated_time_ps < 0.0)
        throw std::invalid_argument("Integrated heat-source time must be finite and nonnegative.");
    const double added = std::accumulate(rates_.begin(), rates_.end(), 0.0) * integrated_time_ps;
    if (!std::isfinite(prescribed_energy_ev() + added))
        throw std::runtime_error("Total heat-source energy ledger overflow.");
    // Validate all cells before updating the ledger.
    for (size_t cell = 0; cell < rates_.size(); ++cell) {
        const double energy = rates_[cell] * integrated_time_ps;
        if (!std::isfinite(pending_[cell] + energy) || !std::isfinite(prescribed_[cell] + energy))
            throw std::runtime_error("Heat-source energy ledger overflow.");
    }
    for (size_t cell = 0; cell < rates_.size(); ++cell) {
        const double energy = rates_[cell] * integrated_time_ps;
        prescribed_[cell] += energy;
        pending_[cell] += energy;
    }
}

double PrescribedHeatSource::deposit_cell(int cell) {
    const double energy = pending_.at(static_cast<size_t>(cell));
    deposited_.at(static_cast<size_t>(cell)) += energy;
    pending_[static_cast<size_t>(cell)] = 0.0;
    return energy;
}

double PrescribedHeatSource::base_power_w() const {
    return std::accumulate(rates_.begin(), rates_.end(), 0.0) / kEvPerWattPs;
}
double PrescribedHeatSource::prescribed_energy_ev() const {
    return std::accumulate(prescribed_.begin(), prescribed_.end(), 0.0);
}
double PrescribedHeatSource::deposited_energy_ev() const {
    return std::accumulate(deposited_.begin(), deposited_.end(), 0.0);
}
double PrescribedHeatSource::pending_energy_ev() const {
    return std::accumulate(pending_.begin(), pending_.end(), 0.0);
}

}  // namespace phonomc
