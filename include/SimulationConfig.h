#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

enum class InitialTemperatureMode { Uniform, Linear };
enum class BoundaryCondition { ThermalReservoir, Periodic, Rough };
enum class HeatSourceProfile { Uniform, Gaussian };
enum class HeatSourceTimeProfile { Constant, Square };
enum class HeatSourceSpectrum { Thermal, Weighted, Gaussian };
enum class InterfaceModel { DiffuseMismatch, AcousticMismatch, MixedMismatch };
// FullMatrix is retained for internal matrix-adapter tests, not a public TOML mode.
enum class CollisionModel { Rta, FullMatrix, Callaway };
inline bool is_linearized_collision(CollisionModel model) {
    return model==CollisionModel::FullMatrix || model==CollisionModel::Callaway;
}

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
const char* to_string(HeatSourceSpectrum spectrum);
const char* to_string(InterfaceModel model);
char boundary_condition_code(BoundaryCondition condition);

struct SimulationConfig {
    std::string input_directory;
    std::string model = "box";
    bool extruded_z = false; // validated constant-height prism; exact clipped grid volumes
    std::vector<double> sizes {200.0, 200.0, 200.0};
    double particle_count = 1e4;
    double time_step = 1.0;
    int iterations = 10000;
    int convergence_write_interval = 10;
    std::uint64_t random_seed = 12345;
    bool compute_kappa = false;
    bool write_cell_heat_flux = false;
    int transport_axis = -1; // -1: legacy grid-based selection; 0/1/2: x/y/z
    std::optional<double> temperature_gradient; // signed K/m, periodic linear-response drive
    // Opt-in gradient-run accelerations; defaults preserve existing trajectories.
    int resample_interval = 0; // steps, 0 disables
    int resample_per_mode_sign = 2; // per cell, mode and deviation sign
    double gradient_warm_start_ps = 0; // homogeneous RTA transient used only as initial guess
    bool convergence_stop = false;
    double convergence_min_time_ps = 2000;
    double convergence_window_ps = 500;
    int convergence_windows = 4;
    double convergence_relative_tolerance = 0.005;
    double convergence_absolute_tolerance = 0.01; // W/(m K)
    bool profile_timers = false;
    bool progress_temperature_summary_only = false;
    bool merge_coplanar_facets = true;  // true: merge coplanar connected triangles into one facet
    double temperature_lookup_dt = 0.1;  // K, lookup-table step for T->E precompute
    double background_temperature = 300.0; // K, invariant deviational-energy reference
    bool lifetime_temperature_is_local = true;
    double lifetime_temperature = 300.0; // K, used when lifetime_temperature_is_local is false
    CollisionModel collision_model = CollisionModel::Rta;
    bool callaway_nonlinear = false; // opt-in displaced Bose targets, conserving finite steps
    bool conservative_boundaries = false; // modal elastic reflection and fixed incoming reservoir flux
    // One converted phono3py energy-generator file per material, in material order.
    std::vector<std::string> scattering_matrix_files;

    GridShape grid;
    InitialTemperatureConfig initial_temperature;
    std::vector<BoundaryCondition> boundary_conditions;
    std::vector<std::string> boundary_position;
    std::vector<double> boundary_values;
    std::vector<std::string> periodic_pair;
    std::string material_folder;
    std::string output_folder;

    // Optional layered multi-material domain.  When material_folders is empty,
    // io.material_folder defines the single legacy material.  Otherwise N
    // folders define N slabs separated by N-1 absolute positions (nm in TOML,
    // converted to Angstrom internally) along material_interface_axis.
    std::vector<std::string> material_folders;
    int material_interface_axis = 0; // 0=x, 1=y, 2=z
    std::vector<double> material_interface_positions;
    InterfaceModel material_interface_model = InterfaceModel::DiffuseMismatch;
    double material_interface_frequency_bin_thz = 0.1;
    double material_interface_amm_fraction = 0.0; // MMM event-level mixing weight
    double material_interface_parallel_bin_inv_a = 0.05;
    std::vector<double> material_mass_densities_kg_m3;


    // Optional local volumetric heat source.
    bool heat_source_enabled = false;
    std::vector<double> heat_source_min;         // uniform region or optional Gaussian clipping box: 3 absolute coordinates in nm input, stored internally as Angstrom
    std::vector<double> heat_source_max;         // uniform region or optional Gaussian clipping box: 3 absolute coordinates in nm input, stored internally as Angstrom
    double heat_source_power_density = 0.0;      // W/m^3 (uniform: region value, gaussian: peak value)
    std::optional<double> heat_source_total_power; // W, mutually exclusive with density inputs
    HeatSourceSpectrum heat_source_spectrum = HeatSourceSpectrum::Thermal;
    std::optional<double> heat_source_frequency_center; // Gaussian modal energy weight, THz
    std::optional<double> heat_source_frequency_sigma;  // standard deviation, THz (not spatial sigma)
    std::vector<int> heat_source_branches; // optional zero-based branch filter
    std::vector<double> heat_source_branch_weights; // weighted spectrum: relative energy per sampled mode
    double heat_source_frequency_min = -1.0;     // optional modal filter, THz for omega/(2*pi); <0 means no lower cutoff
    double heat_source_frequency_max = -1.0;     // optional modal filter, THz for omega/(2*pi); <0 means no upper cutoff
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

    // Optional initial targeted modal over-occupation.
    bool mode_excitation_enabled = false;
    double mode_excitation_frequency_min = -1.0; // THz for omega/(2*pi); <0 means no lower cutoff
    double mode_excitation_frequency_max = -1.0; // THz for omega/(2*pi); <0 means no upper cutoff
    double mode_excitation_occupation_multiplier = 1.0;
};

inline bool uses_linearized_transport(const SimulationConfig& config) {
    return (is_linearized_collision(config.collision_model) && !config.callaway_nonlinear) || config.temperature_gradient.has_value();
}

inline bool uses_conservative_transport(const SimulationConfig& config) {
    return uses_linearized_transport(config) || config.callaway_nonlinear || config.conservative_boundaries;
}

SimulationConfig load_simulation_config(const std::string& path);
void validate_heat_source_config(const SimulationConfig& args);
void validate_scattering_config(const SimulationConfig& args);
std::string create_indexed_output_folder(const std::string& output_folder);
