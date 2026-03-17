#pragma once

#include "SimulationConfig.h"

#include <atomic>
#include <array>
#include <cstdint>
#include <ostream>
#include <random>
#include <utility>
#include <vector>

class SimulationDomain;
class PhononMaterial;

class MonteCarloSolver {
public:
    MonteCarloSolver(const SimulationConfig& args, const SimulationDomain& geometry, const PhononMaterial& phonon);

    void run_timestep();
    int current_timestep() const { return current_timestep_; }
    int openmp_thread_count() const { return openmp_thread_count_; }
    bool profile_timers_enabled() const { return profile_timers_enabled_; }
    void append_profile_summary(std::ostream& out) const;

private:
    using Vec3 = std::array<double, 3>;

    void initialize_particles(const SimulationDomain& geometry, const PhononMaterial& phonon);
    void initialize_particle_modes(const PhononMaterial& phonon);
    void initialize_particle_temperatures(const SimulationDomain& geometry);
    void initialize_particle_velocities(const PhononMaterial& phonon);
    void initialize_reservoir_injection(const SimulationDomain& geometry, const PhononMaterial& phonon);
    void initialize_rough_boundary_scattering(const SimulationDomain& geometry, const PhononMaterial& phonon);
    void initialize_local_heat_source(const SimulationDomain& geometry);
    void apply_local_heat_source();
    std::mt19937_64& thread_rng() const;
    int sample_diffuse_active_mode(int rough_idx, int in_ai, int* source = nullptr) const;
    std::array<int, 2> select_reflected_mode(
        const SimulationDomain& geometry,
        const PhononMaterial& phonon,
        int rough_idx,
        const std::array<int, 2>& in_mode,
        const Vec3& collision_pos,
        double& out_occupation,
        double in_occupation) const;
    void update_collision_cache(const SimulationDomain& geometry, const std::vector<int>& indices);
    void update_collision_cache_single(const SimulationDomain& geometry, int i);
    int nearest_grid_index(const SimulationDomain& geometry, const Vec3& p) const;
    void process_boundary_collision(const SimulationDomain& geometry, int i);
    int remove_absorbed_particles();
    int recover_escaped_particles(const SimulationDomain& geometry);
    std::vector<std::pair<int, double>> inject_particles_from_reservoirs(const SimulationDomain& geometry, const PhononMaterial& phonon);
    void advance_particle(const SimulationDomain& geometry, const PhononMaterial& phonon, int i, double dt_remaining);
    void update_particle_temperatures(const SimulationDomain& geometry, const PhononMaterial& phonon);
    void update_grid_energy_density(const SimulationDomain& geometry, const PhononMaterial& phonon);
    void apply_lifetime_scattering(const PhononMaterial& phonon);
    void update_heat_flux_and_conductivity(const SimulationDomain& geometry);
    void report_timestep_timers_if_needed() const;
    double compute_roughness_specularity(const SimulationDomain& geometry, const PhononMaterial& phonon, int i, int facet) const;
    void write_convergence_header();
    void append_convergence_row() const;
    void write_rough_boundary_mode_map(const SimulationDomain& geometry, const PhononMaterial& phonon) const;
    void ensure_tls_buffers(int thread_count, int nsv);

    static Vec3 add(const Vec3& a, const Vec3& b);
    static Vec3 sub(const Vec3& a, const Vec3& b);
    static Vec3 mul(const Vec3& a, double s);
    static double dot(const Vec3& a, const Vec3& b);
    static double norm(const Vec3& a);
    Vec3 random_unit_vector();

    SimulationConfig args_;
    const SimulationDomain* geometry_ = nullptr;
    const PhononMaterial* phonon_ = nullptr;
    mutable std::mt19937_64 rng_ {std::random_device{}()};
    std::uint64_t rng_seed_base_ = 0;
    int particle_count_ = 0;
    double time_step_ = 1.0;
    double elapsed_time_ = 0.0;
    int current_timestep_ = 0;
    int convergence_write_interval_ = 10;
    const double angstrom_to_meter_ = 1e-10;
    const double evpsa2_to_wm2_ = 1.602176634e13;
    const double wm3_to_evpsa3_ = 6.241509074e-24;
    double particle_density_ = 0.0;
    double push_eps_ = 1e-10;

    std::vector<std::array<int, 2>> particle_modes_;
    std::vector<Vec3> particle_positions_;
    std::vector<Vec3> particle_velocities_;
    std::vector<Vec3> cached_collision_positions_;
    std::vector<double> timesteps_to_collision_;
    std::vector<double> particle_temperatures_;
    std::vector<double> grid_temperatures_;
    std::vector<int> grid_particle_counts_;
    std::vector<double> grid_energy_density_;
    std::vector<double> particle_omega_;
    std::vector<double> particle_occupation_;
    std::vector<double> particle_energies_;
    std::vector<Vec3> grid_heat_flux_;
    double average_heat_flux_along_axis_ = 0.0;
    double thermal_conductivity_ = 0.0;
    double thermal_conductivity_fit_ = 0.0;
    double thermal_conductivity_endpoints_ = 0.0;
    std::vector<int> particle_grid_id_;
    std::vector<int> cached_collision_facets_;
    std::vector<char> cached_collision_conditions_;
    std::vector<std::uint8_t> particle_alive_flags_;

    int reservoir_count_ = 0;
    std::vector<int> reservoir_facets_;
    std::vector<int> facet_to_reservoir_index_;
    std::vector<double> reservoir_temperatures_;
    std::vector<double> reservoir_areas_;
    std::vector<Vec3> reservoir_normals_;
    std::vector<std::array<int, 2>> reservoir_modes_;
    std::vector<std::vector<double>> reservoir_entry_probability_;
    std::vector<int> reservoir_leaving_prev_step_;
    std::vector<int> reservoir_leaving_curr_step_;

    struct RoughFacetData {
        int facet = -1;
        std::vector<double> specularity;
        std::vector<int> spec_match_active;
        std::vector<double> diffuse_creation_rate;
        std::vector<double> diffuse_creation_prob;
        std::vector<int> outgoing_active;
        std::vector<int> outgoing_sorted_active;
        std::vector<double> outgoing_sorted_omega;
        std::vector<int> diffuse_begin;
        std::vector<int> diffuse_end;
        std::vector<int> diffuse_roulette_active;
        std::vector<double> diffuse_roulette_cdf;
    };
    std::vector<int> facet_to_rough_data_;
    std::vector<RoughFacetData> rough_boundary_data_;

    bool local_heat_source_enabled_ = false;
    Vec3 local_heat_source_min_ {0.0, 0.0, 0.0};
    Vec3 local_heat_source_max_ {0.0, 0.0, 0.0};
    Vec3 local_heat_source_center_ {0.0, 0.0, 0.0};
    Vec3 local_heat_source_sigma_ {1.0, 1.0, 1.0};
    std::string local_heat_source_profile_ {"uniform"};
    double local_heat_source_power_density_wm3_ = 0.0;
    std::vector<double> local_heat_source_grid_weights_;

    // Reusable per-thread scratch buffers to avoid per-step allocations in OpenMP paths.
    std::vector<double> energy_tls_buffer_;
    std::vector<int> count_tls_buffer_;
    std::vector<Vec3> flux_tls_buffer_;
    int tls_thread_count_ = 0;
    int tls_nsv_ = 0;

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

    // Rough-boundary selection diagnostics (thread-safe counters).
    mutable std::atomic<long long> rough_events_total_ {0};
    mutable std::atomic<long long> rough_specular_selected_ {0};
    mutable std::atomic<long long> rough_diffuse_selected_ {0};
    mutable std::atomic<long long> rough_fallback_missing_rough_data_ {0};
    mutable std::atomic<long long> rough_fallback_missing_spec_match_ {0};
    mutable std::atomic<long long> rough_fallback_outgoing_pool_ {0};
    mutable std::atomic<long long> rough_fallback_global_random_ {0};

    // Escaped-particle recovery diagnostics.
    long long escaped_recovery_events_ = 0;
    long long escaped_recovered_particles_ = 0;
    int escaped_recovery_check_interval_ = 100;

    // Per-step reservoir bookkeeping diagnostics.
    long long step_absorbed_particles_ = 0;
    long long step_injected_particles_ = 0;
    long long step_recovered_particles_ = 0;
    long long step_net_particles_ = 0;  // injected - absorbed
    long long total_absorbed_particles_ = 0;
    long long total_injected_particles_ = 0;
    long long total_recovered_particles_ = 0;
    long long total_net_particles_ = 0;  // cumulative(injected - absorbed)
};
