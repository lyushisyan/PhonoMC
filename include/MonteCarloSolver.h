#pragma once

#include "SimulationConfig.h"
#include "solver/BackgroundCache.h"
#include "solver/BoundaryLedger.h"
#include "solver/CarrierCollisionOperator.h"
#include "solver/ConvergenceMonitor.h"
#include "solver/EnergyLedger.h"
#include "solver/FullMatrixCollisionOperator.h"
#include "solver/GridState.h"
#include "solver/GridStatistics.h"
#include "solver/InterfaceInteraction.h"
#include "solver/LayeredMaterialLayout.h"
#include "solver/NonlinearCallawayCollision.h"
#include "solver/ParticleStorage.h"
#include "solver/PhononSourceOperator.h"
#include "solver/PrescribedHeatSource.h"
#include "solver/ReservoirInjector.h"
#include "solver/ReservoirSampler.h"
#include "solver/ResultWriter.h"
#include "solver/RoughBoundarySampler.h"
#include "solver/RtaCollisionOperator.h"
#include "solver/TemperatureGradientDrive.h"
#include "solver/TransportOperator.h"

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <random>
#include <utility>
#include <vector>

class SimulationDomain;
class PhononMaterial;

// Public simulation facade. TransportOperator owns event sequencing, numerical
// placement and cache bookkeeping; ReservoirInjector owns thermal-bath insertion.
// The solver coordinates both with sources, collisions, compaction and observations.
// See docs/source/architecture.rst for their shared-state and timestep contracts.
class MonteCarloSolver {
public:
    MonteCarloSolver(const SimulationConfig& args,
                     const SimulationDomain& geometry,
                     const PhononMaterial& phonon);
    MonteCarloSolver(const SimulationConfig& args,
                     const SimulationDomain& geometry,
                     const std::vector<std::reference_wrapper<const PhononMaterial>>& materials);

    // Operators borrow solver-owned services; moving a live solver would invalidate them.
    MonteCarloSolver(const MonteCarloSolver&) = delete;
    MonteCarloSolver& operator=(const MonteCarloSolver&) = delete;
    MonteCarloSolver(MonteCarloSolver&&) = delete;
    MonteCarloSolver& operator=(MonteCarloSolver&&) = delete;

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

    double total_lifetime_energy_residual_ev() const {
        return energy_ledger_.total().lifetime_residual_ev;
    }

    double total_energy_balance_residual_ev() const {
        return energy_ledger_.total().balance_residual_ev;
    }

    double average_heat_flux_along_axis() const { return average_heat_flux_along_axis_; }

    double total_gradient_energy_ev() const { return energy_ledger_.total().drive_energy_ev; }

    long long absorbed_particle_count() const { return boundary_ledger_.total().absorbed; }

    long long injected_particle_count() const { return boundary_ledger_.total().injected; }

    long long interface_event_count() const { return transport_->diagnostics().interface_events; }

    long long interface_transmitted_count() const {
        return transport_->diagnostics().interface_transmitted;
    }

    long long interface_reflected_count() const {
        return transport_->diagnostics().interface_reflected;
    }

    void append_profile_summary(std::ostream& out) const;

private:
    void transport_particles();
    void record_observables_if_due();
    using Vec3 = phonomc::Vec3;

    // Initialization.cpp
    void initialize_particles(const SimulationDomain& geometry);
    void initialize_particle_modes();
    void initialize_particle_temperatures(const SimulationDomain& geometry);
    void initialize_particle_velocities();
    void apply_initial_mode_excitation();
    void initialize_material_layout(const SimulationDomain& geometry);

    // InterfaceInitialization.cpp
    void initialize_interface_mode_banks();

    // RoughBoundaryInitialization.cpp
    void initialize_rough_boundary_scattering(const SimulationDomain& geometry);

    // ReservoirInitialization.cpp
    void initialize_reservoir_injection(const SimulationDomain& geometry);

    // HeatSourceDriving.cpp
    void initialize_local_heat_source(const SimulationDomain& geometry);
    void apply_heat_source(double begin_ps, double end_ps);
    double local_heat_source_integrated_time_factor(double time_begin_ps, double time_end_ps) const;

    // Observables.cpp
    void refresh_thermodynamic_state(const SimulationDomain& geometry);
    double compute_total_thermal_energy_ev(const SimulationDomain& geometry) const;
    void update_heat_flux_and_conductivity(const SimulationDomain& geometry);
    void report_timestep_timers_if_needed() const;
    void sample_convergence();

    // Diagnostics.cpp: capture output snapshots and report console progress.
    void report_progress_if_needed() const;
    void write_convergence_header();
    void append_convergence_row() const;
    phonomc::ConvergenceSnapshot convergence_snapshot() const;
    phonomc::SourceSummarySnapshot source_summary_snapshot() const;
    phonomc::SourceLedgerSnapshot source_ledger_snapshot() const;
    phonomc::RunDiagnosticsSnapshot run_diagnostics_snapshot() const;

    // MonteCarloSolver.cpp: orchestration, material access and shared helpers.
    std::mt19937_64& thread_rng() const;
    void apply_phonon_collisions();
    void apply_temperature_gradient(double time_ps);
    const PhononMaterial& material(int material_id) const;
    const PhononMaterial& particle_material(int particle_index) const;
    int material_index_at(const Vec3& position) const;
    double material_energy_weight(int material_id) const;
    std::vector<double> material_energy_weights() const;
    static Vec3 sub(const Vec3& a, const Vec3& b);
    static double dot(const Vec3& a, const Vec3& b);

    SimulationConfig args_;
    phonomc::ResultWriter result_writer_;
    phonomc::ConvergenceMonitor convergence_monitor_;
    const SimulationDomain* geometry_ = nullptr;
    const PhononMaterial* phonon_ = nullptr;
    std::vector<const PhononMaterial*> materials_;
    mutable std::mt19937_64 rng_;
    mutable std::vector<std::mt19937_64> thread_rngs_;
    int initial_particle_count_ = 0;
    double elapsed_time_ = 0.0;
    int current_timestep_ = 0;
    bool incomplete_timestep_ = false;
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
    double thermal_conductivity_fit_ = 0.0;
    double thermal_conductivity_endpoints_ = 0.0;

    phonomc::ReservoirSampler reservoir_sampler_;
    phonomc::BoundaryLedger boundary_ledger_;
    phonomc::EnergyLedger energy_ledger_;

    phonomc::RoughBoundarySampler rough_sampler_;

    phonomc::LayeredMaterialLayout material_layout_;
    phonomc::InterfaceInteraction interface_interaction_;
    std::unique_ptr<phonomc::TransportOperator> transport_;
    std::unique_ptr<phonomc::ReservoirInjector> reservoir_injector_;

    phonomc::PrescribedHeatSource heat_source_;
    phonomc::PhononSourceOperator source_operator_;
    phonomc::TemperatureGradientDrive gradient_drive_;

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

    // Domain-membership validation; legacy recovery output columns remain zero.
    int escaped_recovery_check_interval_ = 100;
};
