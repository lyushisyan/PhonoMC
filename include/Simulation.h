#pragma once

#include "SimulationConfig.h"

#include <array>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

class Geometry;
class Phonon;

class Simulation {
public:
    Simulation(const SimulationConfig& args, const Geometry& geometry, const Phonon& phonon);

    void run_timestep();
    int current_timestep() const { return current_timestep_; }

private:
    using Vec3 = std::array<double, 3>;

    void initialize_particles(const Geometry& geometry, const Phonon& phonon);
    void initialize_particle_modes(const Phonon& phonon);
    void initialize_particle_temperatures(const Geometry& geometry);
    void initialize_particle_velocities(const Phonon& phonon);
    void initialize_reservoir_injection(const Geometry& geometry, const Phonon& phonon);
    void initialize_rough_boundary_scattering(const Geometry& geometry, const Phonon& phonon);
    void initialize_local_heat_source(const Geometry& geometry);
    void apply_local_heat_source();
    int sample_diffuse_active_mode(int rough_idx, int in_ai) const;
    std::array<int, 2> select_reflected_mode(
        const Geometry& geometry,
        const Phonon& phonon,
        int rough_idx,
        const std::array<int, 2>& in_mode,
        const Vec3& collision_pos,
        double& out_occupation,
        double in_occupation) const;
    void update_collision_cache(const Geometry& geometry, const std::vector<int>& indices);
    int classify_subvolume(const Geometry& geometry, const Vec3& p) const;
    void process_collision(const Geometry& geometry, int i);
    void remove_absorbed_particles();
    std::vector<std::pair<int, double>> inject_particles_from_reservoirs(const Geometry& geometry, const Phonon& phonon);
    void advance_particle(const Geometry& geometry, const Phonon& phonon, int i, double dt_remaining);
    void update_particle_temperatures(const Geometry& geometry, const Phonon& phonon);
    void update_subvolume_energy_density(const Geometry& geometry, const Phonon& phonon);
    void apply_lifetime_scattering(const Phonon& phonon);
    void update_heat_flux_and_conductivity(const Geometry& geometry);
    double compute_roughness_specularity(const Geometry& geometry, const Phonon& phonon, int i, int facet) const;
    void write_convergence_header();
    void append_convergence_row() const;

    static Vec3 add(const Vec3& a, const Vec3& b);
    static Vec3 sub(const Vec3& a, const Vec3& b);
    static Vec3 mul(const Vec3& a, double s);
    static double dot(const Vec3& a, const Vec3& b);
    static double norm(const Vec3& a);
    Vec3 random_unit_vector();

    SimulationConfig args_;
    const Geometry* geometry_ = nullptr;
    const Phonon* phonon_ = nullptr;
    mutable std::mt19937_64 rng_ {std::random_device{}()};
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
    std::vector<double> subvolume_temperatures_;
    std::vector<int> subvolume_particle_counts_;
    std::vector<double> subvolume_energy_density_;
    std::vector<double> particle_omega_;
    std::vector<double> particle_occupation_;
    std::vector<double> particle_energies_;
    std::vector<Vec3> subvolume_heat_flux_;
    double average_heat_flux_along_axis_ = 0.0;
    double thermal_conductivity_ = 0.0;
    double thermal_conductivity_fit_ = 0.0;
    double thermal_conductivity_endpoints_ = 0.0;
    std::vector<int> particle_subvolume_id_;
    std::vector<int> cached_collision_facets_;
    std::vector<char> cached_collision_conditions_;
    std::vector<std::uint8_t> particle_alive_flags_;

    int reservoir_count_ = 0;
    std::vector<int> reservoir_facets_;
    std::vector<double> reservoir_temperatures_;
    std::vector<double> reservoir_areas_;
    std::vector<Vec3> reservoir_normals_;
    std::vector<std::array<int, 2>> reservoir_modes_;
    std::vector<std::vector<double>> reservoir_entry_probability_;
    std::vector<std::vector<double>> reservoir_emit_counter_;

    struct RoughFacetData {
        int facet = -1;
        std::vector<double> specularity;
        std::vector<int> spec_match_active;
        std::vector<int> outgoing_active;
        std::vector<int> outgoing_sorted_active;
        std::vector<double> outgoing_sorted_omega;
        std::vector<int> diffuse_begin;
        std::vector<int> diffuse_end;
    };
    std::vector<int> facet_to_rough_data_;
    std::vector<RoughFacetData> rough_boundary_data_;

    bool local_heat_source_enabled_ = false;
    Vec3 local_heat_source_min_ {0.0, 0.0, 0.0};
    Vec3 local_heat_source_max_ {0.0, 0.0, 0.0};
    double local_heat_source_power_density_wm3_ = 0.0;
    std::vector<std::uint8_t> local_heat_source_subvolume_mask_;
};
