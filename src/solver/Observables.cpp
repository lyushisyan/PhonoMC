#include "MonteCarloSolver.h"
#include "SimulationDomain.h"
#include "solver/StageTimer.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#ifdef PHONOMC_USE_OPENMP
#include <omp.h>
#endif

void MonteCarloSolver::refresh_thermodynamic_state(const SimulationDomain& geometry) {
    phonomc::StageTimer timer(timer_update_temp_, profile_timers_enabled_ && incomplete_timestep_);
    ++temperature_refresh_count_;
    grid_statistics_.refresh_temperatures(particles_,
                                          grid_,
                                          geometry.grid_volumes(),
                                          materials_,
                                          background_cache_,
                                          uses_linearized_transport(args_));
}

double MonteCarloSolver::compute_total_thermal_energy_ev(const SimulationDomain& geometry) const {
    return grid_statistics_.total_thermal_energy(grid_, geometry.grid_volumes(), materials_);
}

void MonteCarloSolver::sample_convergence() {
    if (!args_.convergence_stop) {
        return;
    }
    std::vector<double> values{thermal_conductivity_endpoints_};
    for (const auto& q : grid_.heat_flux) {
        values.push_back(-q[args_.transport_axis] / (*args_.temperature_gradient));
    }
    convergence_monitor_.sample(elapsed_time_, values);
}

void MonteCarloSolver::update_heat_flux_and_conductivity(const SimulationDomain& geometry) {
    grid_statistics_.update_heat_flux(
        particles_, grid_, geometry.grid_volumes(), background_cache_);
    const int nsv = static_cast<int>(grid_.material_ids.size());
    const auto& volumes = geometry.grid_volumes();

    int axis = 0;
    if (args_.transport_axis >= 0) {
        axis = args_.transport_axis;
    } else {
        const int nx = args_.grid.nx;
        const int ny = args_.grid.ny;
        const int nz = args_.grid.nz;
        if (ny > nx && ny >= nz) {
            axis = 1;
        } else if (nz > nx && nz > ny) {
            axis = 2;
        }
    }

    double phi_weighted = 0.0;
    double total_volume = 0.0;
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for reduction(+ : phi_weighted, total_volume)
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const double volume = volumes[static_cast<size_t>(sv)];
        phi_weighted += grid_.heat_flux[static_cast<size_t>(sv)][axis] * volume;
        total_volume += volume;
    }
    average_heat_flux_along_axis_ = phi_weighted / total_volume;

    if (!args_.compute_kappa) {
        thermal_conductivity_fit_ = 0.0;
        thermal_conductivity_endpoints_ = 0.0;
        return;
    }

    if (args_.temperature_gradient) {
        thermal_conductivity_endpoints_ =
            -average_heat_flux_along_axis_ / *args_.temperature_gradient;
        thermal_conductivity_fit_ = 0.0;
        return;
    }
    double T_left = grid_.temperatures.front();
    double T_right = grid_.temperatures.back();
    const auto& rf = geometry.reservoir_facets();
    if (rf.size() >= 2) {
        const auto& fcent = geometry.mesh().facet_centroids();
        int f_left = rf.front();
        int f_right = rf.front();
        for (int f : rf) {
            if (fcent[static_cast<size_t>(f)][axis] < fcent[static_cast<size_t>(f_left)][axis]) {
                f_left = f;
            }
            if (fcent[static_cast<size_t>(f)][axis] > fcent[static_cast<size_t>(f_right)][axis]) {
                f_right = f;
            }
        }
        T_left = geometry.reservoir_value_for_facet(f_left, T_left);
        T_right = geometry.reservoir_value_for_facet(f_right, T_right);
    }

    const auto& bmin = geometry.bounds_min();
    const auto& bmax = geometry.bounds_max();
    const double L = std::abs(bmax[axis] - bmin[axis]) * angstrom_to_meter_;

    // Method 1: endpoint temperature difference gradient.
    const double grad_end = (std::abs(L) > 1e-24) ? ((T_right - T_left) / L) : 0.0;
    if (std::abs(grad_end) > 1e-18) {
        thermal_conductivity_endpoints_ = -average_heat_flux_along_axis_ / grad_end;
    } else {
        thermal_conductivity_endpoints_ = 0.0;
    }

    // Method 2: linear fit over grid centers.
    const auto& centers = geometry.grid_centers();
    double sx = 0.0;
    double st = 0.0;
    double sxx = 0.0;
    double sxt = 0.0;
    int nfit = 0;
    for (int i = 0; i < nsv; ++i) {
        const double x = centers[static_cast<size_t>(i)][axis] * angstrom_to_meter_;
        const double tt = grid_.temperatures[static_cast<size_t>(i)];
        sx += x;
        st += tt;
        sxx += x * x;
        sxt += x * tt;
        ++nfit;
    }
    double grad_fit = 0.0;
    const double dT_window =
        std::abs(grid_.temperatures.back() - grid_.temperatures.front());
    if (nfit >= 2) {
        const double den = static_cast<double>(nfit) * sxx - sx * sx;
        if (std::abs(den) > 1e-30 && dT_window > 1e-6) {
            grad_fit = (static_cast<double>(nfit) * sxt - sx * st) / den;
        }
    }
    if (std::abs(grad_fit) > 1e-18) {
        thermal_conductivity_fit_ = -average_heat_flux_along_axis_ / grad_fit;
    } else {
        thermal_conductivity_fit_ = 0.0;
    }

}
