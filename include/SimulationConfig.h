#pragma once

#include <string>
#include <vector>

struct SimulationConfig {
    std::string input_directory;
    std::string model = "box";
    std::vector<double> sizes {200.0, 200.0, 200.0};
    double particle_count = 1e4;
    double time_step = 1.0;
    int iterations = 10000;
    int convergence_write_interval = 10;
    bool compute_kappa = false;
    bool profile_timers = false;
    bool progress_temperature_summary_only = false;
    double temperature_lookup_dt = 0.1;  // K, lookup-table step for T->E precompute

    std::vector<std::string> grid_layout;
    std::vector<std::string> initial_temperature {"t0", "300"};
    std::vector<std::string> boundary_conditions;
    std::vector<std::string> boundary_position;
    std::vector<double> boundary_values;
    std::vector<std::string> periodic_pair;
    std::string material_folder;
    std::string output_folder;

    // Optional local volumetric heat source.
    bool heat_source_enabled = false;
    std::vector<double> heat_source_min;         // uniform profile: 3 relative values in [0,1]
    std::vector<double> heat_source_max;         // uniform profile: 3 relative values in [0,1]
    double heat_source_power_density = 0.0;      // W/m^3 (uniform: region value, gaussian: peak value)
    std::string heat_source_profile = "uniform"; // uniform | gaussian
    std::vector<double> heat_source_center;      // gaussian profile: 3 relative values in [0,1]
    std::vector<double> heat_source_sigma;       // gaussian profile: 3 relative values in [0,1], <=0 => uniform axis
};

SimulationConfig load_simulation_config(const std::string& path);
std::string create_indexed_output_folder(const std::string& output_folder);
inline SimulationConfig parse_input_file(const std::string& path) { return load_simulation_config(path); }
inline std::string generate_results_folder(const std::string& output_folder) {
    return create_indexed_output_folder(output_folder);
}
