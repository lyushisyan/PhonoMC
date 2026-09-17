#include "MonteCarloSolver.h"
#include "solver/StageTimer.h"
#include "SimulationDomain.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#ifdef PHONOMC_USE_OPENMP
#include <omp.h>
#endif

void MonteCarloSolver::update_particle_temperatures(const SimulationDomain& geometry) {
    phonomc::StageTimer timer(timer_update_temp_, profile_timers_enabled_ && incomplete_timestep_);
    ++temperature_refresh_count_;
    grid_statistics_.refresh_temperatures(particles_, grid_, geometry.grid_volumes(), materials_, background_cache_,
                                          uses_linearized_transport(args_));
}
double MonteCarloSolver::compute_total_thermal_energy_ev(const SimulationDomain& geometry) const {
    return grid_statistics_.total_thermal_energy(grid_, geometry.grid_volumes(), materials_);
}
void MonteCarloSolver::sample_convergence() {
    if (!args_.convergence_stop) return;
    std::vector<double> values{thermal_conductivity_endpoints_};
    for (const auto& q:grid_.heat_flux) values.push_back(-q[args_.transport_axis]/(*args_.temperature_gradient));
    convergence_monitor_.sample(elapsed_time_,values);
}

void MonteCarloSolver::update_heat_flux_and_conductivity(const SimulationDomain& geometry) {
    grid_statistics_.update_heat_flux(particles_, grid_, geometry.grid_volumes(), background_cache_);
    const int nsv = static_cast<int>(grid_.material_ids.size());
    const auto& volumes = geometry.grid_volumes();

    int axis = 0;
    if (args_.transport_axis >= 0) axis = args_.transport_axis;
    else {
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
#pragma omp parallel for reduction(+:phi_weighted,total_volume)
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const double volume = (sv < static_cast<int>(volumes.size()))
            ? std::max(0.0, volumes[static_cast<size_t>(sv)]) : 0.0;
        phi_weighted += grid_.heat_flux[static_cast<size_t>(sv)][axis] * volume;
        total_volume += volume;
    }
    average_heat_flux_along_axis_ = phi_weighted / std::max(1e-30, total_volume);

    if (!args_.compute_kappa) {
        thermal_conductivity_fit_ = 0.0;
        thermal_conductivity_endpoints_ = 0.0;
        thermal_conductivity_ = 0.0;
        return;
    }

    if (args_.temperature_gradient) {
        thermal_conductivity_endpoints_ = -average_heat_flux_along_axis_ / *args_.temperature_gradient;
        thermal_conductivity_ = thermal_conductivity_endpoints_;
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
            if (f < 0 || f >= static_cast<int>(fcent.size())) {
                continue;
            }
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
    int i0 = 0;
    int i1 = nsv - 1;

    double sx = 0.0;
    double st = 0.0;
    double sxx = 0.0;
    double sxt = 0.0;
    int nfit = 0;
    for (int i = i0; i <= i1 && i < static_cast<int>(centers.size()); ++i) {
        const double x = centers[static_cast<size_t>(i)][axis] * angstrom_to_meter_;
        const double tt = grid_.temperatures[static_cast<size_t>(i)];
        sx += x;
        st += tt;
        sxx += x * x;
        sxt += x * tt;
        ++nfit;
    }
    double grad_fit = 0.0;
    const double dT_window = std::abs(
        grid_.temperatures[static_cast<size_t>(std::clamp(i1, 0, nsv - 1))] -
        grid_.temperatures[static_cast<size_t>(std::clamp(i0, 0, nsv - 1))]);
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

    // Backward-compatible scalar uses endpoint gradient method.
    thermal_conductivity_ = thermal_conductivity_endpoints_;
}

// Console-only progress diagnostics; result files are handled by ResultWriter.
void MonteCarloSolver::report_timestep_timers_if_needed() const {
    if (!profile_timers_enabled_ || current_timestep_ <= 0) {
        return;
    }
    if (current_timestep_ % 100 != 0 && current_timestep_ != args_.iterations) {
        return;
    }
    const double denom = std::max(timer_total_, 1e-18);
    auto pct = [denom](double x) { return 100.0 * x / denom; };
    std::cout << "[profile] timesteps=" << current_timestep_
              << " total=" << timer_total_ << "s\n";
    std::cout << "  advance_main=" << pct(timer_advance_main_) << "%\n";
    std::cout << "  remove_absorb_1=" << pct(timer_remove_absorb_1_) << "%\n";
    std::cout << "  inject_build=" << pct(timer_inject_build_) << "%\n";
    std::cout << "  inject_cache=" << pct(timer_inject_cache_) << "%\n";
    std::cout << "  advance_injected=" << pct(timer_advance_injected_) << "%\n";
    std::cout << "  remove_absorb_2=" << pct(timer_remove_absorb_2_) << "%\n";
    std::cout << "  update_temp=" << pct(timer_update_temp_) << "%\n";
    std::cout << "  lifetime=" << pct(timer_lifetime_) << "%\n";
    std::cout << "  stats=" << pct(timer_stats_) << "%\n";
    std::cout << "  source=" << pct(timer_source_) << "%\n";
    std::cout << "  boundary_checks=" << pct(timer_boundary_checks_) << "%\n";
    std::cout << "  energy_ledger=" << pct(timer_energy_ledger_) << "%\n";
}
