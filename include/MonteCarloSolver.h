#pragma once

#include "SimulationConfig.h"
#include "solver/ParticleStorage.h"
#include "solver/GridState.h"
#include "solver/GridStatistics.h"
#include "solver/BackgroundCache.h"
#include "solver/RtaCollisionOperator.h"
#include "solver/CarrierCollisionOperator.h"
#include "solver/NonlinearCallawayCollision.h"
#include "solver/FullMatrixCollisionOperator.h"
#include "solver/TemperatureGradientDrive.h"
#include "solver/PrescribedHeatSource.h"
#include "solver/PhononSourceOperator.h"
#include "solver/ResultWriter.h"
#include "solver/InterfaceModeSampler.h"
#include "solver/AcousticInterfaceSampler.h"
#include "solver/RoughBoundarySampler.h"
#include "solver/ReservoirSampler.h"
#include "solver/BoundaryLedger.h"
#include "solver/EnergyLedger.h"
#include "solver/ConvergenceMonitor.h"

#include <atomic>
#include <array>
#include <cstdint>
#include <functional>
#include <ostream>
#include <random>
#include <utility>
#include <vector>

class SimulationDomain;
class PhononMaterial;

// Public simulation facade. Particle storage and RTA have independent ownership;
// the remaining private methods are implementation partitions in src/solver/.
// See docs/source/architecture.rst for their shared-state and timestep contracts.
class MonteCarloSolver {
public:
    MonteCarloSolver(const SimulationConfig& args, const SimulationDomain& geometry, const PhononMaterial& phonon);
    MonteCarloSolver(
        const SimulationConfig& args,
        const SimulationDomain& geometry,
        const std::vector<std::reference_wrapper<const PhononMaterial>>& materials);

    void run_timestep();
    int current_timestep() const { return current_timestep_; }
    bool converged() const { return convergence_monitor_.converged(); }
    int particle_count() const { return particles_.size(); }
    int openmp_thread_count() const { return openmp_thread_count_; }
    bool profile_timers_enabled() const { return profile_timers_enabled_; }
    double total_heat_source_energy_ev() const { return energy_ledger_.total().source_energy_ev; }
    double prescribed_heat_source_energy_ev() const { return heat_source_.prescribed_energy_ev(); }
    double pending_heat_source_energy_ev() const { return heat_source_.pending_energy_ev(); }
    double total_thermal_energy_ev() const { return energy_ledger_.thermal_energy_ev(); }
    double total_lifetime_energy_residual_ev() const { return energy_ledger_.total().lifetime_residual_ev; }
    double total_energy_balance_residual_ev() const { return energy_ledger_.total().balance_residual_ev; }
    double average_heat_flux_along_axis() const { return average_heat_flux_along_axis_; }
    double total_gradient_energy_ev() const { return energy_ledger_.total().drive_energy_ev; }
    long long absorbed_particle_count() const { return boundary_ledger_.total().absorbed; }
    long long injected_particle_count() const { return boundary_ledger_.total().injected; }
    long long interface_event_count() const { return interface_events_total_.load(std::memory_order_relaxed); }
    long long interface_transmitted_count() const { return interface_transmitted_total_.load(std::memory_order_relaxed); }
    long long interface_reflected_count() const { return interface_reflected_total_.load(std::memory_order_relaxed); }
    void append_profile_summary(std::ostream& out) const;

private:
    using Vec3 = phonomc::Vec3;

    // Initialization.cpp
    void initialize_particles(const SimulationDomain& geometry);
    void initialize_particle_modes();
    void initialize_particle_temperatures(const SimulationDomain& geometry);
    void initialize_particle_velocities();
    void apply_initial_mode_excitation();
    void initialize_material_layout(const SimulationDomain& geometry);
    Vec3 random_unit_vector();

    // Interfaces.cpp
    void initialize_interface_mode_banks();
    double next_interface_time(int particle_index, int& next_material) const;
    void process_material_interface(int particle_index, int next_material);

    // RoughBoundaries.cpp
    void initialize_rough_boundary_scattering(const SimulationDomain& geometry);
    std::array<int, 2> select_reflected_mode(
        const PhononMaterial& phonon,
        int rough_idx,
        const std::array<int, 2>& in_mode,
        double& out_occupation,
        double in_occupation) const;
    double compute_roughness_specularity(const SimulationDomain& geometry, const PhononMaterial& phonon, int i, int facet) const;

    // Transport.cpp
    void update_collision_cache(const SimulationDomain& geometry, const std::vector<int>& indices);
    void update_collision_cache_single(const SimulationDomain& geometry, int i);
    void throw_if_collision_cache_failed(const char* context) const;
    void throw_if_excessive_collisions() const;
    int nearest_grid_index(const SimulationDomain& geometry, const Vec3& p) const;
    void process_boundary_collision(const SimulationDomain& geometry, int i);
    void throw_if_particles_escaped(const SimulationDomain& geometry) const;
    void advance_particle(const SimulationDomain& geometry, int i, double dt_remaining);

    // Reservoirs.cpp
    void initialize_reservoir_injection(const SimulationDomain& geometry);
    std::vector<std::pair<int, double>> inject_particles_from_reservoirs(const SimulationDomain& geometry);

    // HeatSource.cpp
    void initialize_local_heat_source(const SimulationDomain& geometry);
    void apply_local_heat_source_to_occupations(double integrated_time_factor_ps);
    double local_heat_source_integrated_time_factor(double time_begin_ps, double time_end_ps) const;

    // Observables.cpp
    void update_particle_temperatures(const SimulationDomain& geometry);
    double compute_total_thermal_energy_ev(const SimulationDomain& geometry) const;
    void update_heat_flux_and_conductivity(const SimulationDomain& geometry);
    void report_timestep_timers_if_needed() const;
    void sample_convergence();

    // Diagnostics.cpp: capture immutable output values, never format file columns.
    void write_convergence_header();
    void append_convergence_row() const;
    phonomc::ConvergenceSnapshot convergence_snapshot() const;
    phonomc::SourceSummarySnapshot source_summary_snapshot() const;
    phonomc::SourceLedgerSnapshot source_ledger_snapshot() const;
    phonomc::RunDiagnosticsSnapshot run_diagnostics_snapshot() const;

    // MonteCarloSolver.cpp: orchestration, material access and shared helpers.
    std::mt19937_64& thread_rng() const;
    void apply_lifetime_scattering();
    double background_temperature_reference() const;
    const PhononMaterial& material(int material_id) const;
    const PhononMaterial& particle_material(int particle_index) const;
    int material_index_at(const Vec3& position) const;
    double material_energy_weight(int material_id) const;
    static Vec3 add(const Vec3& a, const Vec3& b);
    static Vec3 sub(const Vec3& a, const Vec3& b);
    static Vec3 mul(const Vec3& a, double s);
    static double dot(const Vec3& a, const Vec3& b);
    static double norm(const Vec3& a);

    SimulationConfig args_;
    phonomc::ResultWriter result_writer_;
    phonomc::ConvergenceMonitor convergence_monitor_;
    std::mt19937_64 resampling_rng_;
    long long resampling_removed_ = 0;
    double resampling_residual_ev_ = 0, resampling_seconds_ = 0;
    const SimulationDomain* geometry_ = nullptr;
    const PhononMaterial* phonon_ = nullptr;
    std::vector<const PhononMaterial*> materials_;
    mutable std::mt19937_64 rng_;
    std::uint64_t rng_seed_base_ = 0;
    mutable std::vector<std::mt19937_64> thread_rngs_;
    int initial_particle_count_ = 0;
    double time_step_ = 1.0;
    double elapsed_time_ = 0.0;
    int current_timestep_ = 0;
    bool incomplete_timestep_ = false;
    int convergence_write_interval_ = 10;
    const double angstrom_to_meter_ = 1e-10;
    double particle_density_ = 0.0;
    double particle_spatial_weight_a3_ = 0.0;
    double push_eps_ = 1e-10;
    bool fixed_lifetime_temperature_ = false;
    double background_temperature_ = 300.0;
    double lifetime_temperature_ = 300.0;

    phonomc::GridState grid_;
    std::vector<phonomc::BackgroundCache> background_cache_;
    double average_heat_flux_along_axis_ = 0.0;
    double thermal_conductivity_ = 0.0;
    double thermal_conductivity_fit_ = 0.0;
    double thermal_conductivity_endpoints_ = 0.0;

    phonomc::ReservoirSampler reservoir_sampler_;
    phonomc::BoundaryLedger boundary_ledger_;
    phonomc::EnergyLedger energy_ledger_;

    phonomc::RoughBoundarySampler rough_sampler_;

    phonomc::InterfaceModeSampler interface_sampler_;
    phonomc::AcousticInterfaceSampler acoustic_interface_sampler_;

    phonomc::PrescribedHeatSource heat_source_;
    phonomc::PhononSourceOperator source_operator_;
    phonomc::TemperatureGradientDrive gradient_drive_;
    void apply_temperature_gradient(double time_ps);

    phonomc::ParticleStorage particles_;
    phonomc::RtaCollisionOperator rta_collision_;
    phonomc::CarrierCollisionOperator carrier_collision_;
    phonomc::NonlinearCallawayCollision nonlinear_callaway_;
    phonomc::FullMatrixCollisionOperator full_matrix_collision_;

    phonomc::GridStatistics grid_statistics_;

    bool profile_timers_enabled_ = false;
    int openmp_thread_count_ = 1;
    double timer_total_ = 0.0;
    double timer_advance_main_ = 0.0;
    double timer_remove_absorb_1_ = 0.0;
    double timer_inject_build_ = 0.0;
    double timer_inject_cache_ = 0.0;
    double timer_advance_injected_ = 0.0;
    double timer_remove_absorb_2_ = 0.0;
    double timer_update_temp_ = 0.0;
    double timer_lifetime_ = 0.0;
    double timer_stats_ = 0.0;
    double timer_source_ = 0.0;
    double timer_boundary_checks_ = 0.0;
    double timer_energy_ledger_ = 0.0;
    long long temperature_refresh_count_ = 0;

    // Rough-boundary selection diagnostics (thread-safe counters).
    mutable std::atomic<long long> rough_events_total_ {0};
    mutable std::atomic<long long> rough_specular_selected_ {0};
    mutable std::atomic<long long> rough_diffuse_selected_ {0};
    mutable std::atomic<long long> rough_residual_window_selected_ {0};
    mutable std::atomic<long long> rough_residual_window_fallback_ {0};
    mutable std::atomic<long long> rough_fallback_missing_rough_data_ {0};
    mutable std::atomic<long long> rough_fallback_missing_spec_match_ {0};
    mutable std::atomic<long long> rough_fallback_outgoing_pool_ {0};
    mutable std::atomic<long long> rough_fallback_global_random_ {0};

    // Internal material-interface diagnostics.
    mutable std::atomic<long long> interface_events_total_ {0};
    mutable std::atomic<long long> interface_transmitted_total_ {0};
    mutable std::atomic<long long> interface_reflected_total_ {0};

    // Collision tracing diagnostics. A failed trace is never hidden by
    // changing the physical velocity direction.
    mutable std::atomic<long long> collision_cache_corrections_total_ {0};
    mutable std::atomic<long long> collision_cache_failures_total_ {0};
    mutable std::atomic<int> excessive_collision_particle_ {-1};

    // Domain-membership validation; legacy recovery output columns remain zero.
    int escaped_recovery_check_interval_ = 100;

};
