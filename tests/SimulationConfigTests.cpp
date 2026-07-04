#include "SimulationConfig.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

namespace {
int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        ++failures;
        std::cerr << "FAIL: " << message << '\n';
    }
}

void check_close(double actual, double expected, const std::string& message) {
    check(std::abs(actual - expected) < 1e-12, message);
}

void expect_throw(const std::function<void()>& fn, const std::string& message) {
    try {
        fn();
        check(false, message);
    } catch (const std::runtime_error&) {
    }
}

fs::path make_test_root() {
    const auto tag = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const fs::path root = fs::temp_directory_path() / ("phonomc-config-tests-" + std::to_string(tag));
    fs::create_directories(root);
    return root;
}

fs::path write_file(const fs::path& root, const std::string& name, const std::string& content) {
    const fs::path path = root / name;
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to create test input: " + path.string());
    }
    out << content;
    return path;
}

void test_toml_strong_types(const fs::path& root) {
    const fs::path path = write_file(root, "valid.toml", R"(
[geometry]
model = "box"
sizes = [10, 20, 30]
[simulation]
grid_xyz = [2, 3, 4]
random_seed = 18446744073709551615
initial_temperature = "linear"
background_temperature_mode = "fixed"
background_temperature = 310
lifetime_temperature_mode = "local"
[boundary]
boundary_position = [
  [-0.01, 0, 0, 0.01, 1, 1],
  [0.99, 0, 0, 1.01, 1, 1],
  [0, -0.01, 0, 1, 0.01, 1],
  [0, 0.99, 0, 1, 1.01, 1]
]
boundary_conditions = ["T", "P", "P", "R"]
boundary_values = [300, 0, 0, 1]
periodic_pair = [
  [0.99, 0, 0, 1.01, 1, 1],
  [0, -0.01, 0, 1, 0.01, 1]
]
[heat_source]
profile = "gaussian"
center = [0.5, 0.5, 0.5]
sigma = [0.1, 0.1, 0.1]
power_density = 1e20
)");

    const SimulationConfig cfg = load_simulation_config(path.string());
    check(cfg.grid.nx == 2 && cfg.grid.ny == 3 && cfg.grid.nz == 4, "grid_xyz maps to GridShape");
    check(cfg.random_seed == UINT64_MAX, "random_seed preserves the full unsigned 64-bit range");
    check(cfg.initial_temperature.mode == InitialTemperatureMode::Linear, "linear initial temperature is typed");
    check(cfg.background_temperature_mode == TemperatureReferenceMode::Fixed, "fixed background mode is typed");
    check(cfg.lifetime_temperature_mode == TemperatureReferenceMode::Local, "local lifetime mode is typed");
    check(cfg.boundary_conditions.size() == 4, "all boundary conditions are parsed");
    check(cfg.boundary_conditions[0] == BoundaryCondition::ThermalReservoir, "T boundary is typed");
    check(cfg.boundary_conditions[1] == BoundaryCondition::Periodic, "P boundary is typed");
    check(cfg.boundary_conditions[2] == BoundaryCondition::Periodic, "paired P boundary is typed");
    check(cfg.boundary_conditions[3] == BoundaryCondition::Rough, "R boundary is typed");
    check(cfg.heat_source_profile == HeatSourceProfile::Gaussian, "heat source profile is typed");
    check(cfg.heat_source_enabled, "complete heat source configuration enables the source");
    check_close(cfg.sizes[0], 100.0, "nm sizes convert to Angstrom");
}

void test_uniform_input(const fs::path& root) {
    const fs::path numeric = write_file(root, "numeric.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
initial_temperature = 325.5
)");
    const SimulationConfig numeric_cfg = load_simulation_config(numeric.string());
    check(numeric_cfg.initial_temperature.mode == InitialTemperatureMode::Uniform, "numeric temperature selects uniform mode");
    check_close(numeric_cfg.initial_temperature.uniform_temperature, 325.5, "numeric initial temperature is retained");
}

void test_removed_formats_are_rejected(const fs::path& root) {
    const fs::path legacy = write_file(root, "legacy.txt", "--grid_xyz 4 5 6\n");
    expect_throw([&] { (void) load_simulation_config(legacy.string()); }, "legacy text input is rejected");

    const fs::path top_level = write_file(root, "top-level.toml", R"(
grid_xyz = [1, 1, 1]
)");
    expect_throw([&] { (void) load_simulation_config(top_level.string()); }, "top-level key aliases are rejected");

    const fs::path unknown = write_file(root, "unknown.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
callaway = false
)");
    expect_throw([&] { (void) load_simulation_config(unknown.string()); }, "unknown keys are rejected");

    const fs::path old_temperature = write_file(root, "old-temperature.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
initial_temperature = ["t0", "300"]
)");
    expect_throw([&] { (void) load_simulation_config(old_temperature.string()); }, "old initial-temperature arrays are rejected");

    const fs::path old_boolean = write_file(root, "old-boolean.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
compute_kappa = yes
)");
    expect_throw([&] { (void) load_simulation_config(old_boolean.string()); }, "non-TOML boolean aliases are rejected");

    const fs::path quoted_temperature = write_file(root, "quoted-temperature.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
initial_temperature = "300"
)");
    expect_throw([&] { (void) load_simulation_config(quoted_temperature.string()); }, "quoted numeric temperature is rejected");
}

void test_invalid_values(const fs::path& root) {
    const fs::path bad_boundary = write_file(root, "bad-boundary.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
[boundary]
boundary_conditions = ["X"]
)");
    expect_throw([&] { (void) load_simulation_config(bad_boundary.string()); }, "invalid boundary condition is rejected");

    const fs::path bad_grid = write_file(root, "bad-grid.toml", R"(
[simulation]
grid_xyz = [1, 0, 1]
)");
    expect_throw([&] { (void) load_simulation_config(bad_grid.string()); }, "non-positive grid dimension is rejected");

    const fs::path fractional_grid = write_file(root, "fractional-grid.toml", R"(
[simulation]
grid_xyz = [1, 1.5, 1]
)");
    expect_throw([&] { (void) load_simulation_config(fractional_grid.string()); }, "fractional grid dimension is rejected");

    const fs::path bad_profile = write_file(root, "bad-profile.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
[heat_source]
profile = "triangle"
)");
    expect_throw([&] { (void) load_simulation_config(bad_profile.string()); }, "invalid heat source profile is rejected");

    const fs::path bad_seed = write_file(root, "bad-seed.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
random_seed = -1
)");
    expect_throw([&] { (void) load_simulation_config(bad_seed.string()); }, "negative random seed is rejected");

    const fs::path cylinder = write_file(root, "cylinder.toml", R"(
[geometry]
model = "cylinder"
sizes = [100, 10, 32]
[simulation]
grid_xyz = [1, 1, 1]
)");
    expect_throw([&] { (void) load_simulation_config(cylinder.string()); }, "removed cylinder model is rejected");

    const fs::path boundary_mismatch = write_file(root, "boundary-mismatch.toml", R"(
[simulation]
grid_xyz = [1, 1, 1]
[boundary]
boundary_position = [[-0.01, 0, 0, 0.01, 1, 1]]
boundary_conditions = ["T"]
boundary_values = []
)");
    expect_throw([&] { (void) load_simulation_config(boundary_mismatch.string()); }, "boundary array length mismatch is rejected");
}

void test_repository_examples() {
    const fs::path examples = fs::path(PHONOMC_SOURCE_DIR) / "example";
    const SimulationConfig cross = load_simulation_config((examples / "input_cross_100nm.toml").string());
    const SimulationConfig inplane = load_simulation_config((examples / "input_inplane_x1000nm_z100nm_r1nm.toml").string());
    const SimulationConfig finfet = load_simulation_config((examples / "input_finfet_stl_heat1e20.toml").string());
    check(cross.grid.nx == 10 && cross.compute_kappa, "cross-plane example remains compatible");
    check(inplane.boundary_conditions.size() == 6, "in-plane example remains compatible");
    check(finfet.heat_source_profile == HeatSourceProfile::Gaussian, "FinFET example remains compatible");
}

void test_indexed_output_folders(const fs::path& root) {
    const fs::path base = root / "results" / "run";
    const fs::path first = create_indexed_output_folder(base.string());
    const fs::path second = create_indexed_output_folder(base.string());
    check(first.filename() == "run_0", "first output folder uses index zero");
    check(second.filename() == "run_1", "second output folder increments the index");
    check(fs::is_directory(first) && fs::is_directory(second), "indexed output folders are created");
}
}  // namespace

int main() {
    const fs::path root = make_test_root();
    try {
        test_toml_strong_types(root);
        test_uniform_input(root);
        test_removed_formats_are_rejected(root);
        test_invalid_values(root);
        test_repository_examples();
        test_indexed_output_folders(root);
    } catch (const std::exception& ex) {
        ++failures;
        std::cerr << "Unexpected exception: " << ex.what() << '\n';
    }
    std::error_code ec;
    fs::remove_all(root, ec);
    if (failures != 0) {
        std::cerr << failures << " test assertion(s) failed\n";
        return 1;
    }
    std::cout << "All SimulationConfig tests passed\n";
    return 0;
}
