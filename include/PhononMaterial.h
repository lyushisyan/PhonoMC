#pragma once

#include "material/MaterialData.h"
#include "material/CallawayOperator.h"

#include <array>
#include <random>
#include <vector>

struct SimulationConfig;

class PhononMaterial {
public:
    using Vec3 = std::array<double, 3>;
    using Mode = std::array<int, 2>;

    PhononMaterial(const SimulationConfig& args, int mat_index);
    // File-free construction. Owns its tables; no config, environment or logging.
    explicit PhononMaterial(phonomc::MaterialData data, double temperature_lookup_dt = 0.1,
                            bool include_stationary_modes = false);
    // Explicitly test-only: never selected automatically by the in-memory model.
    static PhononMaterial synthetic_for_testing(double temperature_lookup_dt = 0.1);

    int active_mode_count() const { return active_mode_count_; }
    int qpoint_count() const { return qpoint_count_; }
    int branch_count() const { return branch_count_; }
    const phonomc::MaterialMat3& reciprocal_lattice() const { return reciprocal_; }
    const std::array<int,3>& q_mesh() const { return mesh_; }
    const std::vector<Mode>& active_mode_list() const { return active_mode_list_; }
    Mode active_mode_at(int active_index) const;
    int active_index_for_mode(const Mode& mode) const;
    Mode sample_active_mode(std::mt19937_64& rng) const;
    Vec3 mode_group_velocity(const Mode& mode) const;
    double mode_angular_frequency(const Mode& mode) const;
    double mode_wavevector_norm(const Mode& mode) const;
    Vec3 mode_wavevector(const Mode& mode) const; // angular Cartesian wavevector, 1/Angstrom
    phonomc::CallawayData callaway_data(double reference_temperature) const;
    phonomc::CallawayData linearized_rta_data(double reference_temperature) const;
    double mode_frequency_window(const Mode& mode) const;
    int degenerate_partner_branch(const Mode& mode) const;
    double bose_occupation(double temperature, const Mode& mode) const;
    double mode_heat_capacity(double temperature, const Mode& mode) const; // eV/K
    double linearized_bose_occupation(double temperature, double reference_temperature, const Mode& mode) const;
    double mode_energy(double temperature, const Mode& mode) const;
    double mode_lifetime(double temperature, const Mode& mode) const;
    std::array<double,2> mode_callaway_rates(double temperature, const Mode& mode) const;
    // Exact finite-step RTA weight; a zero linewidth gives zero relaxation.
    double mode_relaxation_weight(double temperature, const Mode& mode, double dt) const;
    double crystal_energy_density(double temperature) const;
    double energy_density_from_temperature(double temperature) const;
    double temperature_from_energy_density(double energy_density) const;
    double normalize_to_energy_density(double x) const;
    Vec3 normalize_to_energy_density(const Vec3& x) const;
    double energy_density_normalization() const;
    double zero_point_energy_density() const { return zero_point_energy_density_; }

private:
    static Vec3 random_unit_vector(std::mt19937_64& rng);
    double mode_gamma(double temperature, const Mode& mode) const;
    struct SyntheticTag {};
    PhononMaterial(SyntheticTag, double temperature_lookup_dt);
    void initialize_material_data(phonomc::MaterialData data);
    void build_fallback_modes();
    int flatten_mode_index(const Mode& mode) const;
    void initialize_temperature_lookup();
    static double lerp(double x0, double x1, double y0, double y1, double x);
    static double interp_linear_clamped(const std::vector<double>& xs, const std::vector<double>& ys, double x);
    void build_degenerate_mode_map();
    void build_active_mode_maps();

    int active_mode_count_ = 1;
    int qpoint_count_ = 0;
    int branch_count_ = 3;
    double hbar_ = 6.582119569e-4;   // eV*ps
    double kb_ = 8.617333262145e-5;  // eV/K
    double unit_cell_volume_ = 1.0;   // A^3
    double zero_point_energy_density_ = 0.0;        // eV/A^3
    double temperature_lookup_dt_ = 0.1;            // K
    bool include_stationary_modes_ = false;
    std::vector<double> mode_angular_frequency_data_;
    std::vector<double> mode_wavevector_norm_data_;
    std::vector<double> mode_frequency_window_data_;
    std::vector<Vec3> mode_group_velocity_data_;
    std::vector<int> degenerate_partner_branch_data_;
    std::vector<Mode> active_mode_list_;
    std::vector<int> flat_to_active_index_;
    std::vector<double> temperature_samples_;
    std::vector<double> gamma_table_;  // (nT, nQ, nB) flattened
    std::vector<double> gamma_normal_, gamma_umklapp_;
    std::vector<double> gamma_extra_; // summed extra resistive linewidths, (nT,nQ,nB)
    std::vector<double> gamma_total_; // gamma + gamma_extra_, used by RTA
    std::vector<Vec3> mode_wavevector_data_, qpoint_fractions_;
    phonomc::MaterialMat3 reciprocal_ {};
    std::array<int,3> mesh_ {};
    std::vector<double> energy_lookup_table_;
    std::vector<double> temperature_lookup_table_;
};
