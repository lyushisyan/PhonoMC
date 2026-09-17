#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace phonomc {
// Value snapshots own their data: formatting needs no solver/material/geometry access.
// Capture only between completed parallel regions; no concurrent solver sampling.
struct InterfaceChannelSnapshot {
    std::array<double,2> incident_acoustic_flux{}, transmitted_flux{}, no_overlap_flux{};
    double geometric_overlap_flux=0, capacity_removed_flux=0, maximum_capacity_ratio=1;
};
struct SourceSummarySnapshot {
    std::string empty_cell_policy = "defer_locally";
    bool source_enabled {};
    double source_base_power {};
    double source_prescribed {};
    double source_pending {};
    double source_deposited {};
    std::string heat_source_profile {};
    std::string heat_source_spectrum {};
    std::optional<double> heat_source_frequency_center {};
    std::optional<double> heat_source_frequency_sigma {};
    std::string heat_source_time_profile {};
    std::optional<double> heat_source_total_power {};
    double heat_source_power_density {};
    double heat_source_frequency_min {};
    double heat_source_frequency_max {};
    std::vector<int> heat_source_branches {};
    std::vector<double> heat_source_branch_weights {};
    double heat_source_time_start {};
    double heat_source_time_end {};
    double heat_source_period {};
    double heat_source_on_duration {};
    double heat_source_duty_cycle {};
    double heat_source_amplitude {};
    long long total_heat_source_injected_particles {};
    double total_heat_source_injected_energy_ev {};
    long long step_heat_source_injected_particles {};
    double step_heat_source_injected_energy_ev {};
};

struct ConvergenceSnapshot {
    bool gradient_driven = false;
    long long gradient_updates = 0;
    double gradient_energy_ev = 0;
    int particle_count {};
    std::vector<double> temperatures {};
    // Optional spatial diagnostics, in W/m^2; indices match grid_centers.csv.
    std::vector<std::array<double, 3>> cell_heat_flux {};
    double source_prescribed {};
    double source_pending {};
    double source_deposited {};
    int current_timestep {};
    double elapsed_time {};
    double average_heat_flux_along_axis {};
    double thermal_conductivity_fit {};
    double thermal_conductivity_endpoints {};
    long long step_absorbed_particles {};
    long long step_injected_particles {};
    long long step_recovered_particles {};
    long long step_net_particles {};
    double step_reservoir_absorbed_energy_ev {};
    double step_reservoir_injected_energy_ev {};
    long long step_heat_source_injected_particles {};
    double step_heat_source_injected_energy_ev {};
    double step_lifetime_energy_residual_ev {};
    double step_energy_balance_residual_ev {};
    double total_thermal_energy_ev {};
};

struct RunDiagnosticsSnapshot {
    bool acceleration_enabled = false, stopped_on_convergence = false, convergence_enabled = false;
    int resample_interval = 0, resample_per_sign = 0, convergence_windows = 0;
    long long resampled_removed = 0;
    double resample_residual_ev = 0, resample_seconds = 0, warm_start_ps = 0, final_time_ps = 0;
    double convergence_min_ps = 0, convergence_window_ps = 0, convergence_rtol = 0, convergence_atol = 0;
    double convergence_last_end_ps = 0, convergence_spread = 0;
    std::vector<std::vector<double>> convergence_block_means;
    bool gradient_driven = false;
    long long gradient_updates = 0;
    double gradient_energy_ev = 0;
    std::string reservoir_generation = "one_to_one";
    bool openmp_enabled {};
    SourceSummarySnapshot source {};
    int particle_count {};
    std::size_t cell_index_rebuilds {};
    std::size_t material_count {};
    std::string material_interface_model {};
    double interface_frequency_bin_thz {};
    double interface_amm_fraction {}, interface_parallel_bin_inv_a {};
    std::vector<double> material_mass_densities_kg_m3;
    std::vector<InterfaceChannelSnapshot> interface_channels;
    std::vector<double> material_carrier_volumes_a3;
    bool fixed_flux_reservoirs {};
    int openmp_thread_count {};
    std::uint64_t rng_seed_base {};
    bool profile_timers_enabled {};
    long long collision_cache_corrections_total {};
    long long collision_cache_failures_total {};
    long long temperature_refresh_count {};
    long long rough_events_total {};
    long long rough_specular_selected {};
    long long rough_diffuse_selected {};
    long long rough_residual_window_selected {};
    long long rough_residual_window_fallback {};
    long long rough_fallback_missing_rough_data {};
    long long rough_fallback_missing_spec_match {};
    long long rough_fallback_outgoing_pool {};
    long long rough_fallback_global_random {};
    long long interface_events_total {};
    long long interface_transmitted_total {};
    long long interface_reflected_total {};
    int escaped_recovery_check_interval {};
    long long total_absorbed_particles {};
    long long total_injected_particles {};
    long long total_recovered_particles {};
    long long total_net_particles {};
    double total_reservoir_absorbed_energy_ev {};
    double total_reservoir_injected_energy_ev {};
    int initial_particle_count {};
    double particle_spatial_weight_a3 {};
    double background_temperature {};
    double total_thermal_energy_ev {};
    double total_lifetime_energy_residual_ev {};
    double total_energy_balance_residual_ev {};
    double timer_total {};
    double timer_advance_main {};
    double timer_remove_absorb_1 {};
    double timer_inject_build {};
    double timer_inject_cache {};
    double timer_advance_injected {};
    double timer_remove_absorb_2 {};
    double timer_update_temp {};
    double timer_lifetime {};
    double timer_stats {};
    double timer_source {};
    double timer_boundary_checks {};
    double timer_energy_ledger {};
};

struct SourceLedgerCell {
    std::array<double, 3> center {};
    int material = 0;
    double rate_evps = 0, prescribed_ev = 0, deposited_ev = 0, pending_ev = 0;
};
struct SourceLedgerSnapshot {
    bool enabled = false;
    std::string spectrum;
    std::vector<SourceLedgerCell> cells;
};
}  // namespace phonomc
