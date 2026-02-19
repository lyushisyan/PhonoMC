#pragma once

#include <string>
#include <vector>

struct SimulationConfig {
    std::string input_directory;
    std::string model = "box";
    std::vector<double> dimensions {200.0, 200.0, 200.0};
    double particle_count = 1e4;
    double time_step = 1.0;
    int iterations = 10000;
    bool compute_thermal_conductivity = false;

    std::vector<std::string> subvolume_layout;
    std::vector<std::string> initial_temperature_profile {"cold"};
    std::vector<std::string> boundary_conditions;
    std::vector<std::string> boundary_positions;
    std::vector<double> boundary_values;
    std::vector<std::string> periodic_pair_positions;
    std::string material_folder;
    std::string results_base_folder;

    // Optional local volumetric heat source.
    bool heat_source_enabled = false;
    std::string heat_source_mode = "relative";   // relative | absolute
    std::vector<double> heat_source_min;         // 3 values
    std::vector<double> heat_source_max;         // 3 values
    double heat_source_power_density = 0.0;      // W/m^3
};

SimulationConfig load_simulation_config(const std::string& path);
std::string create_indexed_results_folder(const std::string& results_base_folder);

// Backward-compatible aliases for older integration points.
using Args = SimulationConfig;
inline SimulationConfig parse_input_file(const std::string& path) { return load_simulation_config(path); }
inline std::string generate_results_folder(const std::string& results_base_folder) {
    return create_indexed_results_folder(results_base_folder);
}
