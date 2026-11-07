#include "BoundarySelectors.h"
#include "SimulationDomain.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <utility>

// Serialization only: preserve file fields, units and write order. Domain setup
// invokes this after geometry and boundary assignment have completed.
// 函数说明：输出几何、边界、网格与输入信息摘要文件。
void SimulationDomain::write_domain_summary(const SimulationConfig& args) const {
    if (args.output_folder.empty()) {
        return;
    }
    namespace fs = std::filesystem;
    fs::create_directories(args.output_folder);
    std::ofstream out(fs::path(args.output_folder) / "summary.txt");
    auto join_strings = [](const std::vector<std::string>& v) -> std::string {
        if (v.empty()) {
            return "[]";
        }
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) {
                oss << ", ";
            }
            oss << v[i];
        }
        oss << "]";
        return oss.str();
    };
    auto join_numbers = [](const std::vector<double>& v) -> std::string {
        if (v.empty()) {
            return "[]";
        }
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) {
                oss << ", ";
            }
            oss << v[i];
        }
        oss << "]";
        return oss.str();
    };
    auto scale_numbers = [](const std::vector<double>& v, double scale) -> std::vector<double> {
        std::vector<double> out = v;
        for (double& x : out) {
            x *= scale;
        }
        return out;
    };
    auto region_count = [](const std::vector<std::string>& v) -> int {
        if (v.empty()) {
            return 0;
        }
        return static_cast<int>(
            phonomc::geometry_detail::parse_boundary_selectors(v, "summary_selector").boxes.size());
    };
    auto boundary_selector_info =
        [](const std::vector<std::string>& v) -> std::pair<std::string, int> {
        if (v.empty()) {
            return {"none", 0};
        }
        const auto parsed =
            phonomc::geometry_detail::parse_boundary_selectors(v, "boundary_position");
        return {"region", static_cast<int>(parsed.boxes.size())};
    };

    int npairs = 0;
    for (int i = 0; i < facet_count_; ++i) {
        if (periodic_pair_[i] > i) {
            ++npairs;
        }
    }

    out << "[input]\n";
    out << "model = " << args.model << '\n';
    out << "sizes_nm = [" << (args.sizes.size() > 0 ? args.sizes[0] / 10.0 : 0.0) << ", "
        << (args.sizes.size() > 1 ? args.sizes[1] / 10.0 : 0.0) << ", "
        << (args.sizes.size() > 2 ? args.sizes[2] / 10.0 : 0.0) << "]\n";
    out << "sizes_angstrom = [" << (args.sizes.size() > 0 ? args.sizes[0] : 0.0) << ", "
        << (args.sizes.size() > 1 ? args.sizes[1] : 0.0) << ", "
        << (args.sizes.size() > 2 ? args.sizes[2] : 0.0) << "]\n";
    out << "particle_count = " << args.particle_count << '\n';
    out << "time_step = " << args.time_step << '\n';
    out << "iterations = " << args.iterations << '\n';
    out << "convergence_write_interval = " << args.convergence_write_interval << '\n';
    out << "random_seed = " << args.random_seed << '\n';
    out << "temperature_lookup_dt = " << args.temperature_lookup_dt << '\n';
    out << "background_temperature = " << args.background_temperature << '\n';
    out << "callaway_formulation = " << (args.callaway_nonlinear ? "nonlinear" : "linearized")
        << '\n';
    out << "conservative_boundaries = " << (uses_conservative_transport(args) ? "true" : "false")
        << '\n';
    if (uses_conservative_transport(args)) {
        out << "collision_model = "
            << (args.collision_model == CollisionModel::Callaway ? "callaway"
                : args.collision_model == CollisionModel::Rta    ? "rta"
                                                                 : "full_matrix")
            << '\n';
        out << "collision_reference_temperature = " << args.background_temperature << '\n';
        out << "reservoir_generation = fixed_incoming_flux\n";
        for (size_t i = 0; i < args.scattering_matrix_files.size(); ++i) {
            out << "scattering_matrix_" << i << " = " << args.scattering_matrix_files[i] << '\n';
        }
    }
    out << "lifetime_temperature = ";
    if (args.lifetime_temperature_is_local) {
        out << "local";
    } else {
        out << args.lifetime_temperature;
    }
    out << '\n';
    if (args.transport_axis >= 0) {
        out << "transport_axis = " << "xyz"[args.transport_axis] << '\n';
    }
    if (args.temperature_gradient) {
        out << "temperature_gradient_K_per_m = " << *args.temperature_gradient << '\n';
        out << "conductivity_definition = -mean_heatflux/imposed_gradient\n";
        out << "temperature_output = T0_plus_periodic_correction\n";
        out << "kappa_int_applicable = false\n";
    }
    out << "compute_kappa = " << (args.compute_kappa ? "true" : "false") << '\n';
    out << "profile_timers = " << (args.profile_timers ? "true" : "false") << '\n';
    out << "progress_temperature_summary_only = "
        << (args.progress_temperature_summary_only ? "true" : "false") << '\n';
    out << "merge_coplanar_facets = " << (args.merge_coplanar_facets ? "true" : "false") << '\n';
    out << "initial_temperature = ";
    if (args.initial_temperature.mode == InitialTemperatureMode::Uniform) {
        out << args.initial_temperature.uniform_temperature;
    } else {
        out << to_string(args.initial_temperature.mode);
    }
    out << '\n';
    out << "grid_xyz = [" << args.grid.nx << ", " << args.grid.ny << ", " << args.grid.nz << "]\n";
    out << "boundary_conditions = ";
    for (size_t i = 0; i < args.boundary_conditions.size(); ++i) {
        if (i > 0) {
            out << " ";
        }
        out << to_string(args.boundary_conditions[i]);
    }
    out << '\n';
    out << "boundary_values = " << join_numbers(args.boundary_values) << '\n';
    const auto [boundary_selector_mode, boundary_selector_count] =
        boundary_selector_info(args.boundary_position);
    out << "boundary_position_selector_mode = " << boundary_selector_mode << '\n';
    out << "boundary_position_selector_count = " << boundary_selector_count << '\n';
    out << "boundary_position_tokens = " << join_strings(args.boundary_position) << '\n';
    out << "periodic_pair_region_count = " << region_count(args.periodic_pair) << '\n';
    out << "periodic_pair_tokens = " << join_strings(args.periodic_pair) << '\n';
    out << "material_folder = " << args.material_folder << '\n';
    out << "material_folders = " << join_strings(args.material_folders) << '\n';
    out << "material_interface_axis = " << "xyz"[args.material_interface_axis] << '\n';
    out << "material_interface_positions_nm = "
        << join_numbers(scale_numbers(args.material_interface_positions, 0.1)) << '\n';
    out << "material_interface_model = " << to_string(args.material_interface_model) << '\n';
    out << "material_interface_amm_fraction = "
        << (args.material_interface_model == InterfaceModel::AcousticMismatch
                ? 1.0
                : args.material_interface_amm_fraction)
        << '\n';
    out << "material_interface_parallel_bin_inv_a = " << args.material_interface_parallel_bin_inv_a
        << '\n';
    out << "output_folder = " << args.output_folder << '\n';
    out << "heat_source_enabled = " << (args.heat_source_enabled ? "true" : "false") << '\n';
    out << "heat_source_profile = " << to_string(args.heat_source_profile) << '\n';
    out << "heat_source_min_nm = " << join_numbers(scale_numbers(args.heat_source_min, 0.1))
        << '\n';
    out << "heat_source_max_nm = " << join_numbers(scale_numbers(args.heat_source_max, 0.1))
        << '\n';
    out << "heat_source_center_nm = " << join_numbers(scale_numbers(args.heat_source_center, 0.1))
        << '\n';
    out << "heat_source_centers_nm = " << join_numbers(scale_numbers(args.heat_source_centers, 0.1))
        << '\n';
    out << "heat_source_sigma_nm = " << join_numbers(scale_numbers(args.heat_source_sigma, 0.1))
        << '\n';
    out << "heat_source_half_width_nm = "
        << join_numbers(scale_numbers(args.heat_source_half_width, 0.1)) << '\n';
    out << "heat_source_min_angstrom = " << join_numbers(args.heat_source_min) << '\n';
    out << "heat_source_max_angstrom = " << join_numbers(args.heat_source_max) << '\n';
    out << "heat_source_center_angstrom = " << join_numbers(args.heat_source_center) << '\n';
    out << "heat_source_centers_angstrom = " << join_numbers(args.heat_source_centers) << '\n';
    out << "heat_source_sigma_angstrom = " << join_numbers(args.heat_source_sigma) << '\n';
    out << "heat_source_half_width_angstrom = " << join_numbers(args.heat_source_half_width)
        << '\n';
    out << "heat_source_power_density = " << args.heat_source_power_density << '\n';
    out << "heat_source_power_densities = " << join_numbers(args.heat_source_power_densities)
        << '\n';
    out << "heat_source_frequency_min_THz = " << args.heat_source_frequency_min << '\n';
    out << "heat_source_frequency_max_THz = " << args.heat_source_frequency_max << '\n';
    out << "heat_source_time_profile = " << to_string(args.heat_source_time_profile) << '\n';
    out << "mode_excitation_enabled = " << (args.mode_excitation_enabled ? "true" : "false")
        << '\n';
    out << "mode_excitation_frequency_min_THz = " << args.mode_excitation_frequency_min << '\n';
    out << "mode_excitation_frequency_max_THz = " << args.mode_excitation_frequency_max << '\n';
    out << "mode_excitation_occupation_multiplier = " << args.mode_excitation_occupation_multiplier
        << '\n';
    out << '\n';

    out << "[geometry]\n";
    out << "volume = " << volume_ << '\n';
    out << "bounds_min = [" << bounds_min_[0] << ", " << bounds_min_[1] << ", " << bounds_min_[2]
        << "]\n";
    out << "bounds_max = [" << bounds_max_[0] << ", " << bounds_max_[1] << ", " << bounds_max_[2]
        << "]\n";
    out << "n_faces = " << mesh_.face_count() << '\n';
    out << "n_facets = " << facet_count_ << '\n';
    out << "n_simplices = " << mesh_.simplex_count() << '\n';
    out << '\n';

    out << "[grid]\n";
    out << "count = " << grid_count_ << '\n';
    out << "lookup = voxel_to_active_grid_o1\n";
    out << '\n';

    out << "[boundary]\n";
    out << "reservoir_facets = " << reservoir_facets_.size() << '\n';
    out << "rough_facets = " << rough_facets_.size() << '\n';
    out << "periodic_pairs = " << npairs << '\n';
    out << '\n';

    out << "[files]\n";
    out << "grid_centers_csv = " << (fs::path(args.output_folder) / "grid_centers.csv").string()
        << '\n';

    std::ofstream centers_out(fs::path(args.output_folder) / "grid_centers.csv");
    centers_out << "index,x_nm,y_nm,z_nm,volume_nm3\n";
    const size_t n = grid_centers_.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& c = grid_centers_[i];
        const double v = (i < grid_volumes_.size()) ? grid_volumes_[i] : 0.0;
        // Internal length unit is Angstrom; export centers in nm for user-facing outputs.
        centers_out << i << "," << (c[0] / 10.0) << "," << (c[1] / 10.0) << "," << (c[2] / 10.0)
                    << "," << (v / 1000.0) << "\n";
    }
}
