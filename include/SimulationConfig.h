#pragma once

#include <cstdint>
#include <string>
#include <vector>

enum class InitialTemperatureMode { Uniform, Linear };
enum class BoundaryCondition { ThermalReservoir, Periodic, Rough };
enum class HeatSourceProfile { Uniform, Gaussian };
enum class HeatSourceTimeProfile { Constant, Square };

struct GridShape {
    int nx = 1;
    int ny = 1;
    int nz = 1;
};

struct InitialTemperatureConfig {
    InitialTemperatureMode mode = InitialTemperatureMode::Uniform;
    double uniform_temperature = 300.0;
};

const char* to_string(InitialTemperatureMode mode);
const char* to_string(BoundaryCondition condition);
const char* to_string(HeatSourceProfile profile);
const char* to_string(HeatSourceTimeProfile profile);
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
    double background_temperature = 300.0; // K, invariant deviational-energy reference
    bool lifetime_temperature_is_local = true;
    double lifetime_temperature = 300.0; // K, used when lifetime_temperature_is_local is false

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
    std::vector<double> heat_source_min;         // uniform region or optional Gaussian clipping box: 3 absolute coordinates in nm input, stored internally as Angstrom
    std::vector<double> heat_source_max;         // uniform region or optional Gaussian clipping box: 3 absolute coordinates in nm input, stored internally as Angstrom
    double heat_source_power_density = 0.0;      // W/m^3 (uniform: region value, gaussian: peak value)
    HeatSourceProfile heat_source_profile = HeatSourceProfile::Uniform;
    std::vector<double> heat_source_center;      // gaussian profile: 3 absolute coordinates in nm input, stored internally as Angstrom
    std::vector<double> heat_source_centers;     // optional gaussian multi-source centers: flattened triples, nm input, stored internally as Angstrom
    std::vector<double> heat_source_sigma;       // gaussian profile: 3 absolute widths in nm input, stored internally as Angstrom; <=0 => uniform axis
    std::vector<double> heat_source_half_width;  // optional local top-hat half widths for sigma<=0 axes, nm input, stored internally as Angstrom
    std::vector<double> heat_source_power_densities; // optional per-source peak values, W/m^3
    HeatSourceTimeProfile heat_source_time_profile = HeatSourceTimeProfile::Constant;
    double heat_source_time_start = 0.0;         // ps
    double heat_source_time_end = -1.0;          // ps, <0 means no end cutoff
    double heat_source_period = 0.0;             // ps, square-wave period
    double heat_source_on_duration = -1.0;       // ps, square-wave on duration; <0 uses duty_cycle
    double heat_source_duty_cycle = 0.5;         // square-wave duty cycle if on_duration is omitted
    double heat_source_amplitude = 1.0;          // multiplier during the on state
};

SimulationConfig load_simulation_config(const std::string& path);
std::string create_indexed_output_folder(const std::string& output_folder);
