#include "ConfigInternal.h"
#include "ConfigText.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
using phonomc::config_detail::to_lower;
using phonomc::config_detail::trim;

// 函数说明：校验温度参考配置，避免非法固定温度进入热物性查询。
void validate_temperature_reference_config(const SimulationConfig& args) {
    if (!std::isfinite(args.background_temperature) || args.background_temperature < 0.0) {
        throw std::runtime_error("background_temperature must be a finite non-negative value.");
    }
    if (!std::isfinite(args.lifetime_temperature) || args.lifetime_temperature < 0.0) {
        throw std::runtime_error("lifetime_temperature must be a finite non-negative value.");
    }
}

size_t flattened_region_count(const std::vector<std::string>& values,
                              const std::string& field_name) {
    if (values.empty()) {
        return 0;
    }
    if (values.front() != "relative_box" || (values.size() - 1) % 6 != 0) {
        throw std::runtime_error(field_name + " has an invalid internal region representation.");
    }
    return (values.size() - 1) / 6;
}

void validate_boundary_config(const SimulationConfig& args) {
    const size_t selector_count =
        flattened_region_count(args.boundary_position, "boundary.boundary_position");
    const size_t periodic_region_count =
        flattened_region_count(args.periodic_pair, "boundary.periodic_pair");
    if (selector_count == 0) {
        if (!args.boundary_conditions.empty() || !args.boundary_values.empty() ||
            periodic_region_count != 0) {
            throw std::runtime_error("boundary.boundary_position is required when boundary "
                                     "conditions, values, or periodic pairs are provided.");
        }
        return;
    }
    if (args.boundary_conditions.size() != selector_count) {
        throw std::runtime_error("boundary.boundary_conditions must contain exactly one entry per "
                                 "boundary_position region.");
    }
    if (args.boundary_values.size() != selector_count) {
        throw std::runtime_error("boundary.boundary_values must contain exactly one entry per "
                                 "boundary_position region; use 0 for periodic boundaries.");
    }

    size_t periodic_condition_count = 0;
    for (size_t i = 0; i < selector_count; ++i) {
        const BoundaryCondition condition = args.boundary_conditions[i];
        const double value = args.boundary_values[i];
        if (condition == BoundaryCondition::ThermalReservoir && value < 0.0) {
            throw std::runtime_error(
                "Thermal reservoir boundary values must be non-negative temperatures.");
        }
        if (condition == BoundaryCondition::Rough && value < 0.0) {
            throw std::runtime_error(
                "Rough boundary values must be non-negative roughness values.");
        }
        if (condition == BoundaryCondition::Periodic) {
            ++periodic_condition_count;
            if (value != 0.0) {
                throw std::runtime_error("Periodic boundary values must be 0.");
            }
        }
    }
    if (periodic_region_count != periodic_condition_count || periodic_region_count % 2 != 0) {
        throw std::runtime_error("boundary.periodic_pair must contain every periodic boundary "
                                 "region exactly once, in pairs.");
    }
}

void validate_mode_excitation_config(const SimulationConfig& args) {
    if (!args.mode_excitation_enabled) {
        return;
    }
    if (!(args.mode_excitation_occupation_multiplier > 0.0) ||
        !std::isfinite(args.mode_excitation_occupation_multiplier)) {
        throw std::runtime_error(
            "mode_excitation.occupation_multiplier must be finite and positive.");
    }
    if (!(args.mode_excitation_frequency_min < 0.0) &&
        !(std::isfinite(args.mode_excitation_frequency_min) &&
          args.mode_excitation_frequency_min >= 0.0)) {
        throw std::runtime_error(
            "mode_excitation.frequency_min must be finite and non-negative if provided.");
    }
    if (!(args.mode_excitation_frequency_max < 0.0) &&
        !(std::isfinite(args.mode_excitation_frequency_max) &&
          args.mode_excitation_frequency_max >= 0.0)) {
        throw std::runtime_error(
            "mode_excitation.frequency_max must be finite and non-negative if provided.");
    }
    if (args.mode_excitation_frequency_min >= 0.0 && args.mode_excitation_frequency_max >= 0.0 &&
        args.mode_excitation_frequency_max < args.mode_excitation_frequency_min) {
        throw std::runtime_error(
            "mode_excitation.frequency_max must be no smaller than frequency_min.");
    }
}

void validate_material_config(const SimulationConfig& args) {
    if (!std::isfinite(args.material_interface_frequency_bin_thz) ||
        args.material_interface_frequency_bin_thz <= 0) {
        throw std::runtime_error(
            "materials.interface_frequency_bin_thz must be finite and positive.");
    }
    if (!std::isfinite(args.material_interface_amm_fraction) ||
        args.material_interface_amm_fraction < 0 || args.material_interface_amm_fraction > 1) {
        throw std::runtime_error("materials.interface_amm_fraction must be in [0,1].");
    }
    if (!std::isfinite(args.material_interface_parallel_bin_inv_a) ||
        args.material_interface_parallel_bin_inv_a <= 0) {
        throw std::runtime_error(
            "materials.interface_parallel_bin_inv_a must be finite and positive.");
    }
    if (args.material_interface_model == InterfaceModel::DiffuseMismatch &&
        args.material_interface_amm_fraction != 0) {
        throw std::runtime_error("DMM cannot have nonzero AMM fraction.");
    }
    if (args.material_interface_model == InterfaceModel::AcousticMismatch ||
        (args.material_interface_model == InterfaceModel::MixedMismatch &&
         args.material_interface_amm_fraction > 0)) {
        if (args.material_folders.size() < 2 ||
            args.material_mass_densities_kg_m3.size() != args.material_folders.size()) {
            throw std::runtime_error("AMM/MMM requires one mass density per material layer.");
        }
        for (double rho : args.material_mass_densities_kg_m3) {
            if (!std::isfinite(rho) || rho <= 0) {
                throw std::runtime_error("Material mass densities must be finite and positive.");
            }
        }
    }
    if (args.material_folders.empty()) {
        return;
    }
    if (args.material_folders.size() < 2) {
        throw std::runtime_error("materials.folders must contain at least two folders; use "
                                 "io.material_folder for a single material.");
    }
    if (args.material_interface_positions.size() + 1 != args.material_folders.size()) {
        throw std::runtime_error(
            "materials.interface_positions must contain exactly folders.size()-1 values.");
    }
    if (!std::is_sorted(args.material_interface_positions.begin(),
                        args.material_interface_positions.end()) ||
        std::adjacent_find(args.material_interface_positions.begin(),
                           args.material_interface_positions.end()) !=
            args.material_interface_positions.end()) {
        throw std::runtime_error("materials.interface_positions must be strictly increasing.");
    }
    const double extent = args.sizes.size() == 3
                              ? args.sizes[static_cast<size_t>(args.material_interface_axis)]
                              : 0.0;
    const int grid_cells = args.material_interface_axis == 0   ? args.grid.nx
                           : args.material_interface_axis == 1 ? args.grid.ny
                                                               : args.grid.nz;
    const double cell_width = extent / static_cast<double>(grid_cells);
    const double alignment_tolerance = std::max(1e-10, 1e-10 * extent);
    for (double position : args.material_interface_positions) {
        if (!std::isfinite(position)) {
            throw std::runtime_error("materials.interface_positions must be finite.");
        }
        // STL bounds can be translated and are only known after loading the mesh.
        // initialize_material_layout checks physical bounds and grid alignment.
        if (to_lower(args.model) != "box") {
            continue;
        }
        if (!(position > 0.0) || !(position < extent)) {
            throw std::runtime_error("materials.interface_positions must lie strictly inside the "
                                     "domain along interface_axis.");
        }
        const double grid_coordinate = position / cell_width;
        if (std::abs(grid_coordinate - std::round(grid_coordinate)) >
            alignment_tolerance / cell_width) {
            throw std::runtime_error(
                "materials.interface_positions must coincide with grid-cell boundaries.");
        }
    }
    for (const std::string& folder : args.material_folders) {
        if (trim(folder).empty()) {
            throw std::runtime_error("materials.folders entries must not be empty.");
        }
    }
}

} // namespace

void validate_heat_source_config(const SimulationConfig& args) {
    if (!args.heat_source_enabled) {
        return;
    }
    if (args.heat_source_total_power.has_value()) {
        if (!std::isfinite(*args.heat_source_total_power) ||
            !(*args.heat_source_total_power > 0.0)) {
            throw std::runtime_error("heat_source.total_power must be finite and positive (W).");
        }
        if (args.heat_source_power_density != 0.0 || !args.heat_source_power_densities.empty()) {
            throw std::runtime_error("heat_source.total_power cannot be combined with "
                                     "power_density or power_densities.");
        }
    } else if (!(args.heat_source_power_density > 0.0) ||
               !std::isfinite(args.heat_source_power_density)) {
        throw std::runtime_error(
            "enabled heat_source requires positive total_power or power_density.");
    }
    for (const auto* values : {&args.heat_source_min,
                               &args.heat_source_max,
                               &args.heat_source_center,
                               &args.heat_source_centers,
                               &args.heat_source_sigma,
                               &args.heat_source_half_width}) {
        for (double value : *values) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("Heat-source spatial parameters must be finite.");
            }
        }
    }
    std::vector<int> branches = args.heat_source_branches;
    std::sort(branches.begin(), branches.end());
    if ((!branches.empty() && branches.front() < 0) ||
        std::adjacent_find(branches.begin(), branches.end()) != branches.end()) {
        throw std::runtime_error("heat_source.branches must contain unique nonnegative indices.");
    }
    if (args.heat_source_spectrum == HeatSourceSpectrum::Gaussian) {
        if (!args.heat_source_frequency_center || !args.heat_source_frequency_sigma ||
            !std::isfinite(*args.heat_source_frequency_center) ||
            *args.heat_source_frequency_center < 0 ||
            !std::isfinite(*args.heat_source_frequency_sigma) ||
            *args.heat_source_frequency_sigma <= 0) {
            throw std::runtime_error("Gaussian spectrum requires finite frequency_center >= 0 and "
                                     "frequency_sigma > 0, in THz.");
        }
    } else if (args.heat_source_frequency_center || args.heat_source_frequency_sigma) {
        throw std::runtime_error(
            "frequency_center/frequency_sigma require heat_source.spectrum='gaussian'.");
    }
    if (args.heat_source_spectrum != HeatSourceSpectrum::Weighted) {
        if (!args.heat_source_branch_weights.empty()) {
            throw std::runtime_error("branch_weights requires heat_source.spectrum='weighted'.");
        }
    } else {
        bool positive = false;
        for (size_t branch = 0; branch < args.heat_source_branch_weights.size(); ++branch) {
            const double weight = args.heat_source_branch_weights[branch];
            if (!std::isfinite(weight) || weight < 0.0) {
                throw std::runtime_error(
                    "heat_source.branch_weights must be finite and nonnegative.");
            }
            if (branches.empty() ||
                std::binary_search(branches.begin(), branches.end(), static_cast<int>(branch))) {
                positive = positive || weight > 0.0;
            }
        }
        if (!positive) {
            throw std::runtime_error(
                "weighted heat_source needs positive weights on selected branches.");
        }
    }
    for (double q : args.heat_source_power_densities) {
        if (!(q > 0.0) || !std::isfinite(q)) {
            throw std::runtime_error(
                "heat_source.power_densities values must be finite and positive.");
        }
    }
    if (!(args.heat_source_frequency_min < 0.0) &&
        !(std::isfinite(args.heat_source_frequency_min) && args.heat_source_frequency_min >= 0.0)) {
        throw std::runtime_error(
            "heat_source.frequency_min must be finite and non-negative if provided.");
    }
    if (!(args.heat_source_frequency_max < 0.0) &&
        !(std::isfinite(args.heat_source_frequency_max) && args.heat_source_frequency_max >= 0.0)) {
        throw std::runtime_error(
            "heat_source.frequency_max must be finite and non-negative if provided.");
    }
    if (args.heat_source_frequency_min >= 0.0 && args.heat_source_frequency_max >= 0.0 &&
        args.heat_source_frequency_max < args.heat_source_frequency_min) {
        throw std::runtime_error(
            "heat_source.frequency_max must be no smaller than frequency_min.");
    }
    if (args.heat_source_profile == HeatSourceProfile::Uniform) {
        if (args.heat_source_min.size() != 3 || args.heat_source_max.size() != 3) {
            throw std::runtime_error(
                "uniform heat_source requires min and max arrays with 3 values.");
        }
    } else {
        const bool has_single_center = args.heat_source_center.size() == 3;
        const bool has_multi_centers = !args.heat_source_centers.empty();
        if ((!has_single_center && !has_multi_centers) || args.heat_source_sigma.size() != 3) {
            throw std::runtime_error("gaussian heat_source requires center=[x,y,z] or "
                                     "centers=[[x,y,z],...] and sigma=[sx,sy,sz].");
        }
        if (has_multi_centers && (args.heat_source_centers.size() % 3 != 0)) {
            throw std::runtime_error("heat_source.centers must contain one or more [x,y,z] rows.");
        }
        const size_t nsrc = has_multi_centers ? args.heat_source_centers.size() / 3 : 1;
        if (!args.heat_source_power_densities.empty() &&
            args.heat_source_power_densities.size() != nsrc) {
            throw std::runtime_error(
                "heat_source.power_densities must match the number of heat_source.centers.");
        }
        if (!args.heat_source_half_width.empty() && args.heat_source_half_width.size() != 3) {
            throw std::runtime_error("heat_source.half_width must contain 3 values if provided.");
        }
        if ((!args.heat_source_min.empty() || !args.heat_source_max.empty()) &&
            (args.heat_source_min.size() != 3 || args.heat_source_max.size() != 3)) {
            throw std::runtime_error("gaussian heat_source min/max clipping box must be omitted or "
                                     "contain 3 values each.");
        }
    }
    if (!std::isfinite(args.heat_source_amplitude) || args.heat_source_amplitude < 0.0) {
        throw std::runtime_error("heat_source.amplitude must be finite and non-negative.");
    }
    if (!std::isfinite(args.heat_source_on_duration)) {
        throw std::runtime_error("heat_source.on_duration must be finite.");
    }
    if (!std::isfinite(args.heat_source_time_start)) {
        throw std::runtime_error("heat_source.start_time must be finite.");
    }
    if (!std::isfinite(args.heat_source_time_end)) {
        throw std::runtime_error("heat_source.end_time must be finite.");
    }
    if (!std::isfinite(args.heat_source_duty_cycle) || args.heat_source_duty_cycle < 0.0 ||
        args.heat_source_duty_cycle > 1.0) {
        throw std::runtime_error("heat_source.duty_cycle must be in [0, 1].");
    }
    if (args.heat_source_time_profile == HeatSourceTimeProfile::Square) {
        if (!(args.heat_source_period > 0.0) || !std::isfinite(args.heat_source_period)) {
            throw std::runtime_error(
                "square heat_source.time_profile requires a finite positive period.");
        }
        if (args.heat_source_on_duration >= 0.0 &&
            (!std::isfinite(args.heat_source_on_duration) ||
             args.heat_source_on_duration > args.heat_source_period)) {
            throw std::runtime_error(
                "heat_source.on_duration must be finite and no larger than period.");
        }
    }
}


void validate_scattering_config(const SimulationConfig& args) {
    if (args.callaway_nonlinear &&
        (args.collision_model != CollisionModel::Callaway || args.temperature_gradient)) {
        throw std::runtime_error(
            "Nonlinear Callaway requires model=callaway and finite-domain thermal driving");
    }
    if (args.resample_interval != 0) {
        throw std::runtime_error(
            "Solver resampling is unsupported: sampling.resample_interval must be 0. "
            "Existing-carrier merging changes quadrature weights.");
    }
    if (args.resample_per_mode_sign < 1 ||
        args.resample_per_mode_sign > 1024 || !std::isfinite(args.gradient_warm_start_ps) ||
        args.gradient_warm_start_ps < 0) {
        throw std::runtime_error("Invalid gradient sampling or warm-start settings.");
    }
    if ((args.gradient_warm_start_ps > 0 || args.convergence_stop) &&
        !args.temperature_gradient) {
        throw std::runtime_error(
            "Sampling, warm start and convergence stopping currently require gradient driving.");
    }
    if (!std::isfinite(args.convergence_min_time_ps) || args.convergence_min_time_ps < 0 ||
        !std::isfinite(args.convergence_window_ps) || args.convergence_window_ps <= 0 ||
        args.convergence_windows < 4 || args.convergence_windows > 100 ||
        !std::isfinite(args.convergence_relative_tolerance) ||
        args.convergence_relative_tolerance <= 0 ||
        !std::isfinite(args.convergence_absolute_tolerance) ||
        args.convergence_absolute_tolerance < 0) {
        throw std::runtime_error(
            "Invalid convergence windows or tolerances (at least four windows required).");
    }
    if (args.convergence_stop &&
        (!args.compute_kappa ||
         args.convergence_window_ps < 4 * args.time_step * args.convergence_write_interval)) {
        throw std::runtime_error("Convergence stopping requires compute_kappa and at least four "
                                 "output intervals per window.");
    }
    if (args.transport_axis < -1 || args.transport_axis > 2) {
        throw std::runtime_error("Invalid transport_axis.");
    }
    if (args.temperature_gradient) {
        if (!std::isfinite(*args.temperature_gradient) || *args.temperature_gradient == 0 ||
            args.transport_axis < 0) {
            throw std::runtime_error("Temperature-gradient driving requires a finite nonzero "
                                     "gradient and explicit simulation.transport_axis.");
        }
        if (args.model != "box" || args.material_folders.size() > 1 ||
            !args.material_interface_positions.empty()) {
            throw std::runtime_error(
                "Temperature-gradient driving currently requires a single-material box.");
        }
        if (args.collision_model == CollisionModel::FullMatrix) {
            throw std::runtime_error("Temperature-gradient driving supports RTA or Callaway.");
        }
        if (!(args.background_temperature > 0) || !std::isfinite(args.background_temperature) ||
            args.lifetime_temperature_is_local ||
            args.lifetime_temperature != args.background_temperature) {
            throw std::runtime_error("Temperature-gradient driving requires lifetime_temperature = "
                                     "background_temperature > 0.");
        }
        if (args.initial_temperature.mode != InitialTemperatureMode::Uniform ||
            args.initial_temperature.uniform_temperature != args.background_temperature ||
            args.mode_excitation_enabled || args.heat_source_enabled) {
            throw std::runtime_error(
                "Temperature-gradient driving requires uniform initial_temperature = "
                "background_temperature and no heat source or mode excitation.");
        }
        for (auto condition : args.boundary_conditions) {
            if (condition == BoundaryCondition::ThermalReservoir) {
                throw std::runtime_error(
                    "Temperature-gradient driving cannot use thermal reservoirs.");
            }
        }
    }
    if (args.collision_model == CollisionModel::Callaway) {
        if (!args.scattering_matrix_files.empty()) {
            throw std::runtime_error("Callaway uses gamma_N/gamma_U, not scattering.matrix_files.");
        }
        if (!std::isfinite(args.background_temperature) || args.background_temperature <= 0) {
            throw std::runtime_error("Callaway requires a positive fixed background_temperature.");
        }
        return;
    }
    if (args.collision_model == CollisionModel::Rta) {
        if (!args.scattering_matrix_files.empty()) {
            throw std::runtime_error(
                "scattering.matrix_files is not supported by the public rta/callaway models.");
        }
        return;
    }
    if (args.collision_model != CollisionModel::FullMatrix) {
        throw std::runtime_error("Invalid collision model.");
    }
    if (!std::isfinite(args.background_temperature) || args.background_temperature <= 0.0) {
        throw std::runtime_error(
            "Full-matrix scattering requires a positive fixed background_temperature.");
    }
    const size_t count = args.material_folders.empty() ? 1 : args.material_folders.size();
    if (args.scattering_matrix_files.size() != count) {
        throw std::runtime_error("scattering.matrix_files must contain one matrix per material.");
    }
    for (const auto& file : args.scattering_matrix_files) {
        if (file.empty()) {
            throw std::runtime_error("Scattering matrix paths must not be empty.");
        }
    }
}

void phonomc::config_detail::validate_loaded_config(const SimulationConfig& args) {
    if (!(args.temperature_lookup_dt > 0.0) || !std::isfinite(args.temperature_lookup_dt)) {
        throw std::runtime_error("temperature_lookup_dt must be a finite positive value.");
    }
    validate_scattering_config(args);
    validate_boundary_config(args);
    validate_heat_source_config(args);
    validate_mode_excitation_config(args);
    validate_material_config(args);
    validate_simulation_controls(args);
    if (args.grid.nx <= 0 || args.grid.ny <= 0 || args.grid.nz <= 0) {
        throw std::runtime_error("grid_xyz values must be positive integers.");
    }
    if (args.model == "box") {
        if (args.sizes.size() != 3 ||
            std::any_of(args.sizes.begin(), args.sizes.end(), [](double value) {
                return !std::isfinite(value) || !(value > 0.0);
            })) {
            throw std::runtime_error("geometry.sizes for box must contain three positive values.");
        }
    }
}

void validate_simulation_controls(const SimulationConfig& args) {
    validate_temperature_reference_config(args);
    if (args.initial_temperature.mode == InitialTemperatureMode::Uniform &&
        (!std::isfinite(args.initial_temperature.uniform_temperature) ||
         args.initial_temperature.uniform_temperature < 0.0)) {
        throw std::runtime_error("initial_temperature must be a finite non-negative temperature.");
    }
    if (args.convergence_write_interval <= 0) {
        throw std::runtime_error("convergence_write_interval must be a positive integer.");
    }
    // ParticleStorage uses int indices. Round only after the representable
    // population is established; never silently replace an empty population.
    const double count = std::round(args.particle_count);
    if (!std::isfinite(count) || count < 1 || count > std::numeric_limits<int>::max()) {
        throw std::runtime_error("particle_count must round to an integer in [1, INT_MAX].");
    }
    if (!(args.time_step > 0.0) || !std::isfinite(args.time_step)) {
        throw std::runtime_error("time_step must be a finite positive value.");
    }
    if (args.iterations <= 0) {
        throw std::runtime_error("iterations must be a positive integer.");
    }
}
