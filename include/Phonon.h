#pragma once

#include "SimulationConfig.h"

#include <array>
#include <string>
#include <random>
#include <vector>

class Phonon {
public:
    using Vec3 = std::array<double, 3>;
    using Mode = std::array<int, 2>;

    Phonon(const SimulationConfig& args, int mat_index);

    int active_mode_count() const { return active_mode_count_; }
    int qpoint_count() const { return qpoint_count_; }
    int branch_count() const { return branch_count_; }
    const std::vector<Mode>& active_mode_list() const { return active_mode_list_; }
    Mode active_mode_at(int active_index) const;
    int active_index_for_mode(const Mode& mode) const;
    Mode sample_active_mode(std::mt19937_64& rng) const;
    Vec3 mode_group_velocity(const Mode& mode) const;
    double mode_angular_frequency(const Mode& mode) const;
    double mode_wavevector_norm(const Mode& mode) const;
    double mode_frequency_window(const Mode& mode) const;
    int degenerate_partner_branch(const Mode& mode) const;
    double bose_occupation(double temperature, const Mode& mode) const;
    double mode_energy(double temperature, const Mode& mode) const;
    double mode_lifetime(double temperature, const Mode& mode) const;
    double crystal_energy_density(double temperature) const;
    double temperature_from_energy_density(double energy_density) const;
    double normalize_to_energy_density(double x) const;
    Vec3 normalize_to_energy_density(const Vec3& x) const;

private:
    static Vec3 random_unit_vector(std::mt19937_64& rng);
    static int nearest_index(const std::vector<double>& arr, double x);
    bool load_hdf5_data(const SimulationConfig& args, int mat_index, std::string* err);
    void build_fallback_modes(const SimulationConfig& args);
    int flatten_mode_index(const Mode& mode) const;
    void initialize_temperature_lookup();
    bool load_poscar_lattice_volume(const std::string& folder, std::string* err);
    static double lerp(double x0, double x1, double y0, double y1, double x);
    static double interp_linear_clamped(const std::vector<double>& xs, const std::vector<double>& ys, double x);
    void build_degenerate_mode_map();

    int active_mode_count_ = 1;
    int qpoint_count_ = 0;
    int branch_count_ = 3;
    double hbar_ = 6.582119569e-4;   // eV*ps
    double kb_ = 8.617333262145e-5;  // eV/K
    double unit_cell_volume_ = 1.0;   // A^3
    double zero_point_energy_density_ = 0.0;        // eV/A^3
    std::array<std::array<double, 3>, 3> reciprocal_lattice_ {{
        std::array<double, 3>{1.0, 0.0, 0.0},
        std::array<double, 3>{0.0, 1.0, 0.0},
        std::array<double, 3>{0.0, 0.0, 1.0}
    }};
    std::array<int, 3> mesh_grid_size_ {1, 1, 1};
    std::string material_folder_path_;
    std::vector<Vec3> qpoint_fractions_;
    std::vector<double> mode_angular_frequency_data_;
    std::vector<double> mode_wavevector_norm_data_;
    std::vector<double> mode_frequency_window_data_;
    std::vector<Vec3> mode_group_velocity_data_;
    std::vector<int> degenerate_partner_branch_data_;
    std::vector<Mode> active_mode_list_;
    std::vector<int> flat_to_active_index_;
    std::vector<int> active_to_flat_index_;
    std::vector<double> temperature_samples_;
    std::vector<double> gamma_table_;  // (nT, nQ, nB) flattened
    std::vector<double> energy_lookup_table_;
    std::vector<double> temperature_lookup_table_;
};
