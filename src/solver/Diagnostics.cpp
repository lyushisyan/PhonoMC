#include "MonteCarloSolver.h"
#include "SimulationDomain.h"
#include <stdexcept>

// Snapshot collection is the only output adapter allowed to read solver state.

phonomc::SourceSummarySnapshot MonteCarloSolver::source_summary_snapshot() const {
    phonomc::SourceSummarySnapshot s;
    if (args_.collision_model == CollisionModel::FullMatrix) s.empty_cell_policy="create_modal_carriers";
    s.source_enabled = heat_source_.enabled();
    s.source_base_power = heat_source_.base_power_w();
    s.source_prescribed = heat_source_.prescribed_energy_ev();
    s.source_pending = heat_source_.pending_energy_ev();
    s.source_deposited = heat_source_.deposited_energy_ev();
    s.heat_source_profile = to_string(args_.heat_source_profile);
    s.heat_source_spectrum = to_string(args_.heat_source_spectrum);
    s.heat_source_frequency_center = args_.heat_source_frequency_center;
    s.heat_source_frequency_sigma = args_.heat_source_frequency_sigma;
    s.heat_source_time_profile = to_string(args_.heat_source_time_profile);
    s.heat_source_total_power = args_.heat_source_total_power;
    s.heat_source_power_density = args_.heat_source_power_density;
    s.heat_source_frequency_min = args_.heat_source_frequency_min;
    s.heat_source_frequency_max = args_.heat_source_frequency_max;
    s.heat_source_branches = args_.heat_source_branches;
    s.heat_source_branch_weights = args_.heat_source_branch_weights;
    s.heat_source_time_start = args_.heat_source_time_start;
    s.heat_source_time_end = args_.heat_source_time_end;
    s.heat_source_period = args_.heat_source_period;
    s.heat_source_on_duration = args_.heat_source_on_duration;
    s.heat_source_duty_cycle = args_.heat_source_duty_cycle;
    s.heat_source_amplitude = args_.heat_source_amplitude;
    s.total_heat_source_injected_particles = energy_ledger_.total().source_updates;
    s.total_heat_source_injected_energy_ev = energy_ledger_.total().source_energy_ev;
    s.step_heat_source_injected_particles = energy_ledger_.step().source_updates;
    s.step_heat_source_injected_energy_ev = energy_ledger_.step().source_energy_ev;
    return s;
}

phonomc::ConvergenceSnapshot MonteCarloSolver::convergence_snapshot() const {
    phonomc::ConvergenceSnapshot s;
    s.gradient_driven = args_.temperature_gradient.has_value();
    s.gradient_updates = energy_ledger_.step().drive_updates;
    s.gradient_energy_ev = energy_ledger_.step().drive_energy_ev;
    s.particle_count = particles_.size();
    s.temperatures = grid_.temperatures;
    if (s.gradient_driven || args_.write_cell_heat_flux) {
        s.cell_heat_flux.reserve(grid_.heat_flux.size());
        for (const auto& q : grid_.heat_flux) s.cell_heat_flux.push_back({q[0], q[1], q[2]});
    }
    s.source_prescribed = heat_source_.prescribed_energy_ev();
    s.source_pending = heat_source_.pending_energy_ev();
    s.source_deposited = heat_source_.deposited_energy_ev();
    s.current_timestep = current_timestep_;
    s.elapsed_time = elapsed_time_;
    s.average_heat_flux_along_axis = average_heat_flux_along_axis_;
    s.thermal_conductivity_fit = thermal_conductivity_fit_;
    s.thermal_conductivity_endpoints = thermal_conductivity_endpoints_;
    s.step_absorbed_particles = boundary_ledger_.step().absorbed;
    s.step_injected_particles = boundary_ledger_.step().injected;
    s.step_recovered_particles = 0;
    s.step_net_particles = boundary_ledger_.step().net();
    s.step_reservoir_absorbed_energy_ev = boundary_ledger_.step().absorbed_energy_ev;
    s.step_reservoir_injected_energy_ev = boundary_ledger_.step().injected_energy_ev;
    s.step_heat_source_injected_particles = energy_ledger_.step().source_updates;
    s.step_heat_source_injected_energy_ev = energy_ledger_.step().source_energy_ev;
    s.step_lifetime_energy_residual_ev = energy_ledger_.step().lifetime_residual_ev;
    s.step_energy_balance_residual_ev = energy_ledger_.step().balance_residual_ev;
    s.total_thermal_energy_ev = energy_ledger_.thermal_energy_ev();
    return s;
}

phonomc::RunDiagnosticsSnapshot MonteCarloSolver::run_diagnostics_snapshot() const {
    phonomc::RunDiagnosticsSnapshot s;
    s.acceleration_enabled=args_.resample_interval>0 || args_.gradient_warm_start_ps>0 || args_.convergence_stop;
    s.resample_interval=args_.resample_interval;s.resample_per_sign=args_.resample_per_mode_sign;
    s.resampled_removed=resampling_removed_;s.resample_residual_ev=resampling_residual_ev_;
    s.resample_seconds=resampling_seconds_;s.warm_start_ps=args_.gradient_warm_start_ps;
    s.final_time_ps=elapsed_time_;s.stopped_on_convergence=converged();
    s.convergence_enabled=args_.convergence_stop;s.convergence_windows=args_.convergence_windows;
    s.convergence_min_ps=args_.convergence_min_time_ps;s.convergence_window_ps=args_.convergence_window_ps;
    s.convergence_rtol=args_.convergence_relative_tolerance;s.convergence_atol=args_.convergence_absolute_tolerance;
    s.convergence_last_end_ps=convergence_monitor_.last_window_end();s.convergence_spread=convergence_monitor_.maximum_spread();
    for (const auto& b:convergence_monitor_.blocks()) s.convergence_block_means.push_back(b);
    s.gradient_driven = args_.temperature_gradient.has_value();
    s.gradient_updates = energy_ledger_.total().drive_updates;
    s.gradient_energy_ev = energy_ledger_.total().drive_energy_ev;
    if (uses_conservative_transport(args_)) s.reservoir_generation="fixed_incoming_flux";
    s.particle_count = particles_.size();
    s.cell_index_rebuilds = grid_.cells.rebuild_count();
    s.material_count = materials_.size();
    s.material_interface_model = to_string(args_.material_interface_model);
    s.interface_frequency_bin_thz = args_.material_interface_frequency_bin_thz;
    s.interface_amm_fraction = args_.material_interface_model==InterfaceModel::AcousticMismatch ? 1 : args_.material_interface_model==InterfaceModel::MixedMismatch ? args_.material_interface_amm_fraction : 0;
    s.interface_parallel_bin_inv_a = args_.material_interface_parallel_bin_inv_a;
    s.material_mass_densities_kg_m3 = args_.material_mass_densities_kg_m3;
    if(s.interface_amm_fraction>0)for(size_t i=0;i+1<materials_.size();++i){
        const auto& a=acoustic_interface_sampler_.audit(i);
        s.interface_channels.push_back({a.incident_acoustic_flux,a.transmitted_flux,a.no_overlap_flux,
            a.geometric_overlap_flux,a.capacity_removed_flux,a.maximum_capacity_ratio});
    }
    s.material_carrier_volumes_a3 = particles_.material_spatial_volumes_a3;
    s.fixed_flux_reservoirs = uses_conservative_transport(args_) || materials_.size()>1;
    s.openmp_thread_count = openmp_thread_count_;
    s.rng_seed_base = rng_seed_base_;
    s.profile_timers_enabled = profile_timers_enabled_;
    s.collision_cache_corrections_total = collision_cache_corrections_total_.load(std::memory_order_relaxed);
    s.collision_cache_failures_total = collision_cache_failures_total_.load(std::memory_order_relaxed);
    s.temperature_refresh_count = temperature_refresh_count_;
    s.rough_events_total = rough_events_total_.load(std::memory_order_relaxed);
    s.rough_specular_selected = rough_specular_selected_.load(std::memory_order_relaxed);
    s.rough_diffuse_selected = rough_diffuse_selected_.load(std::memory_order_relaxed);
    s.rough_residual_window_selected = rough_residual_window_selected_.load(std::memory_order_relaxed);
    s.rough_residual_window_fallback = rough_residual_window_fallback_.load(std::memory_order_relaxed);
    s.rough_fallback_missing_rough_data = rough_fallback_missing_rough_data_.load(std::memory_order_relaxed);
    s.rough_fallback_missing_spec_match = rough_fallback_missing_spec_match_.load(std::memory_order_relaxed);
    s.rough_fallback_outgoing_pool = rough_fallback_outgoing_pool_.load(std::memory_order_relaxed);
    s.rough_fallback_global_random = rough_fallback_global_random_.load(std::memory_order_relaxed);
    s.interface_events_total = interface_events_total_.load(std::memory_order_relaxed);
    s.interface_transmitted_total = interface_transmitted_total_.load(std::memory_order_relaxed);
    s.interface_reflected_total = interface_reflected_total_.load(std::memory_order_relaxed);
    s.escaped_recovery_check_interval = escaped_recovery_check_interval_;
    s.total_absorbed_particles = boundary_ledger_.total().absorbed;
    s.total_injected_particles = boundary_ledger_.total().injected;
    s.total_recovered_particles = 0;
    s.total_net_particles = boundary_ledger_.total().net();
    s.total_reservoir_absorbed_energy_ev = boundary_ledger_.total().absorbed_energy_ev;
    s.total_reservoir_injected_energy_ev = boundary_ledger_.total().injected_energy_ev;
    s.initial_particle_count = initial_particle_count_;
    s.particle_spatial_weight_a3 = particle_spatial_weight_a3_;
    s.background_temperature = background_temperature_;
    s.total_thermal_energy_ev = energy_ledger_.thermal_energy_ev();
    s.total_lifetime_energy_residual_ev = energy_ledger_.total().lifetime_residual_ev;
    s.total_energy_balance_residual_ev = energy_ledger_.total().balance_residual_ev;
    s.timer_total = timer_total_;
    s.timer_advance_main = timer_advance_main_;
    s.timer_remove_absorb_1 = timer_remove_absorb_1_;
    s.timer_inject_build = timer_inject_build_;
    s.timer_inject_cache = timer_inject_cache_;
    s.timer_advance_injected = timer_advance_injected_;
    s.timer_remove_absorb_2 = timer_remove_absorb_2_;
    s.timer_update_temp = timer_update_temp_;
    s.timer_lifetime = timer_lifetime_;
    s.timer_stats = timer_stats_;
    s.timer_source = timer_source_;
    s.timer_boundary_checks = timer_boundary_checks_;
    s.timer_energy_ledger = timer_energy_ledger_;
#ifdef PHONOMC_USE_OPENMP
    s.openmp_enabled = true;
#endif
    s.source = source_summary_snapshot();
    return s;
}

phonomc::SourceLedgerSnapshot MonteCarloSolver::source_ledger_snapshot() const {
    phonomc::SourceLedgerSnapshot s;
    s.enabled = heat_source_.enabled();
    s.spectrum = to_string(args_.heat_source_spectrum);
    if (!s.enabled) return s;
    const auto& centers = geometry_->grid_centers();
    s.cells.reserve(heat_source_.cell_pending_ev().size());
    for (std::size_t cell = 0; cell < heat_source_.cell_pending_ev().size(); ++cell) {
        s.cells.push_back({centers.at(cell), grid_.material_ids.at(cell),
            heat_source_.cell_rates_evps().at(cell), heat_source_.cell_prescribed_ev().at(cell),
            heat_source_.cell_deposited_ev().at(cell), heat_source_.cell_pending_ev().at(cell)});
    }
    return s;
}
void MonteCarloSolver::write_convergence_header() {
    result_writer_.write_convergence_header(grid_.temperatures.size(),args_.temperature_gradient.has_value(),args_.write_cell_heat_flux);
}
void MonteCarloSolver::append_convergence_row() const {
    if (!result_writer_.enabled()) return;
    result_writer_.append_convergence(convergence_snapshot());
    if (heat_source_.enabled()) result_writer_.write_source_ledger(source_ledger_snapshot());
}
void MonteCarloSolver::append_profile_summary(std::ostream& out) const {
    phonomc::ResultWriter::format_run_summary(out, run_diagnostics_snapshot());
    if (!out) throw std::runtime_error("Failed writing solver summary.");
    if (result_writer_.enabled() && heat_source_.enabled())
        result_writer_.write_source_ledger(source_ledger_snapshot());
}
