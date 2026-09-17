#include "solver/ResultWriter.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace phonomc {
namespace {
template<class Formatter>
void write_file(const std::filesystem::path& path, std::ios::openmode mode, Formatter format) {
    std::ofstream out(path, mode);
    if (!out) throw std::runtime_error("Cannot open output file: " + path.string());
    format(out);
    out.close(); // Detect failures deferred until flushing/closing, not just open().
    if (!out) throw std::runtime_error("Failed writing output file: " + path.string());
}
}  // namespace

void ResultWriter::write_convergence_header(std::size_t cells, bool gradient_driven, bool write_cell_heat_flux) const {
    if (!enabled()) return;
    write_file(std::filesystem::path(folder_) / "convergence.txt", std::ios::trunc,
        [cells, gradient_driven](std::ostream& out) { format_convergence_header(out, cells, gradient_driven); });
    if (gradient_driven || write_cell_heat_flux)
        write_file(std::filesystem::path(folder_) / "grid_heat_flux.csv", std::ios::trunc,
            [](std::ostream& out) {
                out << "timestep,time_ps,cell_index,temperature_K,qx_W_m2,qy_W_m2,qz_W_m2\n";
            });
}
void ResultWriter::append_convergence(const ConvergenceSnapshot& s) const {
    if (!enabled()) return;
    if (!s.cell_heat_flux.empty() && s.cell_heat_flux.size() != s.temperatures.size())
        throw std::invalid_argument("Cell heat-flux and temperature snapshot sizes differ.");
    write_file(std::filesystem::path(folder_) / "convergence.txt", std::ios::app,
        [&](std::ostream& out) { format_convergence(out, s); });
    if (!s.cell_heat_flux.empty())
        write_file(std::filesystem::path(folder_) / "grid_heat_flux.csv", std::ios::app,
            [&](std::ostream& out) {
                out << std::setprecision(17);
                for (std::size_t cell = 0; cell < s.cell_heat_flux.size(); ++cell) {
                    const auto& q = s.cell_heat_flux[cell];
                    out << s.current_timestep << ',' << s.elapsed_time << ',' << cell << ','
                        << s.temperatures[cell] << ',' << q[0] << ',' << q[1] << ',' << q[2] << '\n';
                }
            });
}
void ResultWriter::write_source_ledger(const SourceLedgerSnapshot& s) const {
    if (!enabled() || !s.enabled) return;
    write_file(std::filesystem::path(folder_) / "heat_source_ledger.txt", std::ios::trunc,
        [&](std::ostream& out) { format_source_ledger(out, s); });
}

void ResultWriter::format_convergence_header(std::ostream& out, std::size_t cells, bool gradient_driven) {
    out << "# timestep time_ps";
    for (std::size_t i = 0; i < cells; ++i) {
        out << " T_" << (i + 1);
    }
    out << " heatflux kappa_int kappa_eff particle_count absorbed injected recovered net"
           " reservoir_absorbed_energy_ev reservoir_injected_energy_ev"
           " hs_occupation_updates hs_injected_energy_ev lifetime_energy_residual_ev"
           " energy_balance_residual_ev total_thermal_energy_ev"
           " hs_prescribed_total_ev hs_pending_total_ev hs_ledger_residual_ev";
    if (gradient_driven) out << " gradient_occupation_updates gradient_energy_ev";
    out << '\n';
}

void ResultWriter::format_convergence(std::ostream& out, const ConvergenceSnapshot& s) {
    out << s.current_timestep << " " << s.elapsed_time;
    for (double tsv : s.temperatures) {
        out << " " << tsv;
    }
    out << std::setprecision(12);
    out << " " << s.average_heat_flux_along_axis << " " << s.thermal_conductivity_fit << " " << s.thermal_conductivity_endpoints
        << " " << s.particle_count
        << " " << s.step_absorbed_particles << " " << s.step_injected_particles << " " << s.step_recovered_particles << " " << s.step_net_particles
        << " " << s.step_reservoir_absorbed_energy_ev << " " << s.step_reservoir_injected_energy_ev
        << " " << s.step_heat_source_injected_particles << " " << s.step_heat_source_injected_energy_ev
        << " " << s.step_lifetime_energy_residual_ev << " " << s.step_energy_balance_residual_ev
        << " " << s.total_thermal_energy_ev
        << " " << s.source_prescribed
        << " " << s.source_pending
        << " " << (s.source_prescribed - s.source_deposited -
                   s.source_pending);
    if (s.gradient_driven) out << " " << s.gradient_updates << " " << s.gradient_energy_ev;
    out << '\n';
}

void ResultWriter::format_source_summary(std::ostream& out, const SourceSummarySnapshot& s) {
    out << "\n[local_heat_source_occupation_injection]\n";
    out << "enabled = " << (s.source_enabled ? "true" : "false") << '\n';
    out << "profile = " << s.heat_source_profile << '\n';
    if (s.heat_source_total_power)
        out << "specified_total_power_w = " << *s.heat_source_total_power << '\n';
    else
        out << "power_density_wm3 = " << s.heat_source_power_density << '\n';
    out << "model = prescribed_lattice_power\n";
    out << "electron_transport = false\n";
    out << "empty_cell_policy = " << s.empty_cell_policy << '\n';
    out << "spectrum = " << s.heat_source_spectrum << '\n';
    if (s.heat_source_spectrum == "gaussian") {
        out << "frequency_center_thz = " << s.heat_source_frequency_center.value_or(0) << '\n';
        out << "frequency_sigma_thz = " << s.heat_source_frequency_sigma.value_or(0) << '\n';
        out << "spectral_normalization = existing_carrier_weights_within_each_cell\n";
        out << "spectral_weight = exp(-0.5*((f-f_center)/sigma_f)^2); energy_per_mode\n";
    }
    out << "power_input = " << (s.heat_source_total_power ? "total_power" : "power_density") << '\n';
    out << "base_total_power_w = " << s.source_base_power << '\n';
    out << "frequency_min_thz = " << s.heat_source_frequency_min << '\n';
    out << "frequency_max_thz = " << s.heat_source_frequency_max << '\n';
    out << "branches_zero_based = [";
    for (size_t i = 0; i < s.heat_source_branches.size(); ++i)
        out << (i ? ", " : "") << s.heat_source_branches[i];
    out << "]\nbranch_weights_per_sampled_mode = [";
    for (size_t i = 0; i < s.heat_source_branch_weights.size(); ++i)
        out << (i ? ", " : "") << s.heat_source_branch_weights[i];
    out << "]\n";
    out << "prescribed_energy_total_ev = " << s.source_prescribed << '\n';
    out << "pending_energy_total_ev = " << s.source_pending << '\n';
    out << "deposition_complete = " << (s.source_pending == 0.0 ? "true" : "false") << '\n';
    out << "source_ledger_residual_ev = " << s.source_prescribed -
        s.source_deposited - s.source_pending << '\n';
    out << "time_profile = " << s.heat_source_time_profile << '\n';
    out << "time_start_ps = " << s.heat_source_time_start << '\n';
    out << "time_end_ps = " << s.heat_source_time_end << '\n';
    out << "period_ps = " << s.heat_source_period << '\n';
    out << "on_duration_ps = " << s.heat_source_on_duration << '\n';
    out << "duty_cycle = " << s.heat_source_duty_cycle << '\n';
    out << "amplitude = " << s.heat_source_amplitude << '\n';
    out << "occupation_updates_total = " << s.total_heat_source_injected_particles << '\n';
    out << "injected_energy_total_ev = " << s.total_heat_source_injected_energy_ev << '\n';
    out << "occupation_updates_last_step = " << s.step_heat_source_injected_particles << '\n';
    out << "injected_energy_last_step_ev = " << s.step_heat_source_injected_energy_ev << '\n';
}

void ResultWriter::format_run_summary(std::ostream& out, const RunDiagnosticsSnapshot& s) {
    if (s.acceleration_enabled) {
        out << std::setprecision(17) << "\n[gradient_acceleration]\n"
            << "resample_interval = " << s.resample_interval << '\n'
            << "per_mode_sign = " << s.resample_per_sign << '\n'
            << "removed_carriers_total = " << s.resampled_removed << '\n'
            << "resampling_energy_residual_ev = " << s.resample_residual_ev << '\n'
            << "resampling_seconds = " << s.resample_seconds << '\n'
            << "warm_start_rta_age_ps = " << s.warm_start_ps << '\n'
            << "final_time_ps = " << s.final_time_ps << '\n'
            << "termination = " << (s.stopped_on_convergence?"stationarity_screen_passed":"iteration_limit_or_manual_snapshot") << '\n'
            << "convergence_enabled = " << s.convergence_enabled << '\n'
            << "convergence_min_time_ps = " << s.convergence_min_ps << '\n'
            << "convergence_window_ps = " << s.convergence_window_ps << '\n'
            << "convergence_required_windows = " << s.convergence_windows << '\n'
            << "convergence_relative_tolerance = " << s.convergence_rtol << '\n'
            << "convergence_absolute_tolerance = " << s.convergence_atol << '\n'
            << "convergence_last_window_end_ps = " << s.convergence_last_end_ps << '\n'
            << "convergence_max_spread_W_mK = " << s.convergence_spread << '\n';
        for (size_t i=0;i<s.convergence_block_means.size();++i) {
            out << "window_" << i << "_global_then_cells_W_mK =";
            for (double v:s.convergence_block_means[i]) out << ' ' << v;
            out << '\n';
        }
    }
    if (s.gradient_driven)
        out << "\n[temperature_gradient_accounting]\noccupation_updates = " << s.gradient_updates
            << "\nnet_energy_ev = " << s.gradient_energy_ev << "\n";
    out << std::setprecision(17);
    out << "\n[performance]\n";
    out << "openmp_enabled = " << (s.openmp_enabled ? "true" : "false") << '\n';
    out << "openmp_thread_count = " << s.openmp_thread_count << '\n';
    out << "random_seed = " << s.rng_seed_base << '\n';
    out << "profile_timers_enabled = " << (s.profile_timers_enabled ? "true" : "false") << '\n';
    out << "collision_cache_corrections_total = "
        << s.collision_cache_corrections_total << '\n';
    out << "collision_cache_failures_total = "
        << s.collision_cache_failures_total << '\n';
    out << "transport_failure_policy = strict\n";
    out << "excessive_collision_recoveries_total = 0\n";
    out << "cell_index_rebuilds = " << s.cell_index_rebuilds << '\n';
    out << "temperature_refreshes = " << s.temperature_refresh_count << '\n';
    const long long rough_total = s.rough_events_total;
    const long long rough_spec = s.rough_specular_selected;
    const long long rough_diff = s.rough_diffuse_selected;
    const long long rough_residual_window = s.rough_residual_window_selected;
    const long long rough_residual_fallback = s.rough_residual_window_fallback;
    const long long rough_fb_no_data = s.rough_fallback_missing_rough_data;
    const long long rough_fb_no_spec = s.rough_fallback_missing_spec_match;
    const long long rough_fb_pool = s.rough_fallback_outgoing_pool;
    const long long rough_fb_rand = s.rough_fallback_global_random;
    const double rough_den = std::max(1.0, static_cast<double>(rough_total));
    out << "\n[rough_boundary_diagnostics]\n";
    out << "events_total = " << rough_total << '\n';
    out << "specular_selected = " << rough_spec << '\n';
    out << "diffuse_selected = " << rough_diff << '\n';
    out << "specular_ratio = " << (static_cast<double>(rough_spec) / rough_den) << '\n';
    out << "diffuse_ratio = " << (static_cast<double>(rough_diff) / rough_den) << '\n';
    out << "residual_flux_frequency_window_selected = " << rough_residual_window << '\n';
    out << "residual_flux_frequency_window_fallback = " << rough_residual_fallback << '\n';
    out << "fallback_missing_rough_data = " << rough_fb_no_data << '\n';
    out << "fallback_missing_spec_match = " << rough_fb_no_spec << '\n';
    out << "fallback_outgoing_pool = " << rough_fb_pool << '\n';
    out << "fallback_global_random = " << rough_fb_rand << '\n';
    out << "fallback_total = "
        << (rough_residual_fallback + rough_fb_no_data + rough_fb_no_spec + rough_fb_pool + rough_fb_rand)
        << '\n';
    out << "mode_map_csv = disabled\n";
    const long long interface_events = s.interface_events_total;
    const long long interface_transmitted = s.interface_transmitted_total;
    const long long interface_reflected = s.interface_reflected_total;
    const double interface_den = std::max(1.0, static_cast<double>(interface_events));
    out << "\n[material_interface_diagnostics]\n";
    out << "material_count = " << s.material_count << '\n';
    out << "model = " << s.material_interface_model << '\n';
    if (s.material_count>1) {
        out << "frequency_bin_thz = " << s.interface_frequency_bin_thz << '\n';
        out << "amm_fraction = " << s.interface_amm_fraction << '\n';
        if(s.interface_amm_fraction>0) {
            out << "amm_approximation = dispersive_scalar_acoustic_no_mode_conversion\n";
            out << "amm_optical_treatment = specular_reflection\n";
            out << "amm_kernel = reciprocal_linearized_qcell_overlap_v2\n";
            out << "amm_bin_role = candidate_search_only\n";
            out << "amm_channel_flux_units = inverse_Angstrom2_inverse_ps\n";
            for(size_t i=0;i<s.interface_channels.size();++i){
                const auto& a=s.interface_channels[i];const std::string prefix="interface_"+std::to_string(i)+"_";
                out << prefix << "geometric_overlap_flux = " << a.geometric_overlap_flux << '\n';
                out << prefix << "capacity_removed_flux = " << a.capacity_removed_flux << '\n';
                out << prefix << "maximum_capacity_ratio = " << a.maximum_capacity_ratio << '\n';
                for(int side=0;side<2;++side){const std::string key=prefix+"side_"+std::to_string(side)+"_";
                    out << key << "incident_acoustic_flux = " << a.incident_acoustic_flux[side] << '\n';
                    out << key << "transmitted_flux = " << a.transmitted_flux[side] << '\n';
                    out << key << "no_overlap_flux = " << a.no_overlap_flux[side] << '\n';
                }
            }
            out << "parallel_wavevector_bin_inv_a = " << s.interface_parallel_bin_inv_a << '\n';
            for(size_t m=0;m<s.material_mass_densities_kg_m3.size();++m)
                out << "material_" << m << "_mass_density_kg_m3 = " << s.material_mass_densities_kg_m3[m] << '\n';
        }
        out << "carrier_sampling = equal_physical_state_weight\n";
        out << "fixed_flux_reservoirs = " << (s.fixed_flux_reservoirs ? "true" : "false") << '\n';
        for (size_t mid=0;mid<s.material_carrier_volumes_a3.size();++mid)
            out << "material_" << mid << "_carrier_volume_a3 = " << s.material_carrier_volumes_a3[mid] << '\n';
    }
    out << "events_total = " << interface_events << '\n';
    out << "transmitted_total = " << interface_transmitted << '\n';
    out << "reflected_total = " << interface_reflected << '\n';
    out << "transmission_ratio = "
        << (static_cast<double>(interface_transmitted) / interface_den) << '\n';
    out << "reflection_ratio = "
        << (static_cast<double>(interface_reflected) / interface_den) << '\n';
    out << "\n[escaped_particle_recovery]\n";
    out << "check_interval = " << s.escaped_recovery_check_interval << '\n';
    out << "policy = strict\n";
    out << "recovery_events = 0\n";
    out << "recovered_particles_total = 0\n";
    out << "\n[reservoir_balance_diagnostics]\n";
    out << "reservoir_gen = " << s.reservoir_generation << '\n';
    out << "absorbed_particles_total = " << s.total_absorbed_particles << '\n';
    out << "injected_particles_total = " << s.total_injected_particles << '\n';
    out << "recovered_particles_total = " << s.total_recovered_particles << '\n';
    out << "net_particles_total = " << s.total_net_particles << '\n';
    out << "absorbed_deviational_energy_total_ev = " << s.total_reservoir_absorbed_energy_ev << '\n';
    out << "injected_deviational_energy_total_ev = " << s.total_reservoir_injected_energy_ev << '\n';
    out << "\n[energy_weighting]\n";
    out << "initial_particle_count = " << s.initial_particle_count << '\n';
    out << "current_particle_count = " << s.particle_count << '\n';
    out << "particle_spatial_weight_a3 = " << s.particle_spatial_weight_a3 << '\n';
    out << "background_reference_temperature_k = " << s.background_temperature << '\n';
    out << "total_thermal_energy_ev = " << s.total_thermal_energy_ev << '\n';
    out << "lifetime_energy_residual_total_ev = " << s.total_lifetime_energy_residual_ev << '\n';
    out << "energy_balance_residual_total_ev = " << s.total_energy_balance_residual_ev << '\n';
    format_source_summary(out, s.source);
    if (!s.profile_timers_enabled) {
        return;
    }
    out << "profile_total_seconds = " << s.timer_total << '\n';
    const double denom = std::max(s.timer_total, 1e-18);
    auto write_seg = [&out, denom](const char* name, double sec) {
        out << name << "_seconds = " << sec << '\n';
        out << name << "_percent = " << (100.0 * sec / denom) << '\n';
    };
    write_seg("advance_main", s.timer_advance_main);
    write_seg("remove_absorb_1", s.timer_remove_absorb_1);
    write_seg("inject_build", s.timer_inject_build);
    write_seg("inject_cache", s.timer_inject_cache);
    write_seg("advance_injected", s.timer_advance_injected);
    write_seg("remove_absorb_2", s.timer_remove_absorb_2);
    write_seg("update_temp", s.timer_update_temp);
    write_seg("lifetime", s.timer_lifetime);
    write_seg("stats", s.timer_stats);
    write_seg("source", s.timer_source);
    write_seg("boundary_checks", s.timer_boundary_checks);
    write_seg("energy_ledger", s.timer_energy_ledger);
    const double accounted = s.timer_advance_main + s.timer_remove_absorb_1 +
        s.timer_inject_build + s.timer_inject_cache + s.timer_advance_injected +
        s.timer_remove_absorb_2 + s.timer_update_temp + s.timer_lifetime + s.timer_stats +
        s.timer_source + s.timer_boundary_checks + s.timer_energy_ledger;
    write_seg("other", std::max(0.0, s.timer_total - accounted));
}

void ResultWriter::format_source_ledger(std::ostream& out, const SourceLedgerSnapshot& s) {
    out << std::setprecision(17)
        << "# prescribed_lattice_power; not an electron-phonon calculation\n"
        << "# spectrum=" << s.spectrum
        << "; pending energy is not yet represented by phonons\n"
        << "# cell x_A y_A z_A material base_power_W prescribed_eV deposited_eV pending_eV\n";
    for (std::size_t cell = 0; cell < s.cells.size(); ++cell) {
        const auto& c = s.cells[cell];
        out << cell << " " << c.center[0] << " " << c.center[1] << " " << c.center[2]
            << " " << c.material << " " << c.rate_evps / 6.241509074e6
            << " " << c.prescribed_ev << " " << c.deposited_ev << " " << c.pending_ev << '\n';
    }
}
}  // namespace phonomc
