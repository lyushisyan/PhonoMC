#pragma once

#include <cstdint>
#include <string>
#include <vector>

enum class TemperatureReferenceMode { Local, Fixed };
enum class InitialTemperatureMode { Uniform, Linear };
enum class BoundaryCondition { ThermalReservoir, Periodic, Rough };
enum class HeatSourceProfile { Uniform, Gaussian };

struct GridShape {
    int nx = 1;
    int ny = 1;
    int nz = 1;
};

struct InitialTemperatureConfig {
    InitialTemperatureMode mode = InitialTemperatureMode::Uniform;
    double uniform_temperature = 300.0;
};

const char* to_string(TemperatureReferenceMode mode);
const char* to_string(InitialTemperatureMode mode);
const char* to_string(BoundaryCondition condition);
const char* to_string(HeatSourceProfile profile);
char boundary_condition_code(BoundaryCondition condition);

struct SimulationConfig {
    std::string input_directory;
    std::string model = "box";
    std::vector<double> sizes {200.0, 200.0, 200.0};
    double particle_count = 1e4;
    double time_step = 1.0;
    int iterations = 10000;
    int convergence_write_interval = 10;
    std::uint64_t random_seed = 12345;
    bool compute_kappa = false;
    bool profile_timers = false;
    bool progress_temperature_summary_only = false;
    bool merge_coplanar_facets = true;  // true: merge coplanar connected triangles into one facet
    double temperature_lookup_dt = 0.1;  // K, lookup-table step for T->E precompute
    TemperatureReferenceMode background_temperature_mode = TemperatureReferenceMode::Local;
    double background_temperature = 300.0; // K, used in Fixed mode
    TemperatureReferenceMode lifetime_temperature_mode = TemperatureReferenceMode::Local;
    double lifetime_temperature = 300.0; // K, used in Fixed mode

    GridShape grid;
    InitialTemperatureConfig initial_temperature;
    std::vector<BoundaryCondition> boundary_conditions;
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
    HeatSourceProfile heat_source_profile = HeatSourceProfile::Uniform;
    std::vector<double> heat_source_center;      // gaussian profile: 3 relative values in [0,1]
    std::vector<double> heat_source_sigma;       // gaussian profile: 3 relative values in [0,1], <=0 => uniform axis
};

SimulationConfig load_simulation_config(const std::string& path);
std::string create_indexed_output_folder(const std::string& output_folder);
