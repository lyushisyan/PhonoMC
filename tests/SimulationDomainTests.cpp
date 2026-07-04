#include "SimulationConfig.h"
#include "SimulationDomain.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <random>
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

void expect_throw(const std::function<void()>& fn, const std::string& message) {
    try {
        fn();
        check(false, message);
    } catch (const std::runtime_error&) {
    }
}

fs::path make_test_root() {
    const auto tag = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const fs::path root = fs::temp_directory_path() / ("phonomc-domain-tests-" + std::to_string(tag));
    fs::create_directories(root);
    return root;
}

SimulationConfig cross_config(const fs::path& output) {
    SimulationConfig cfg = load_simulation_config(
        (fs::path(PHONOMC_SOURCE_DIR) / "example" / "input_cross_100nm.toml").string());
    cfg.output_folder = output.string();
    return cfg;
}
}  // namespace

int main() {
    const fs::path root = make_test_root();
    try {
        SimulationConfig valid = cross_config(root / "valid");
        const SimulationDomain domain(valid);
        check(domain.mesh().facet_count() == 6, "box has six merged facets");
        const double expected_box_volume = valid.sizes[0] * valid.sizes[1] * valid.sizes[2];
        check(std::abs(domain.volume() - expected_box_volume) <= expected_box_volume * 1e-12,
              "signed surface volume matches the analytical box volume");
        check(domain.fast_grid_index_enabled(), "box enables constant-time grid lookup");
        for (int i = 0; i < domain.grid_count(); ++i) {
            check(domain.fast_grid_index(domain.grid_centers()[static_cast<size_t>(i)]) == i,
                  "box grid center maps to its exact grid index");
        }
        check(domain.reservoir_facets().size() == 2, "cross-plane input has two reservoirs");
        for (int facet = 0; facet < domain.mesh().facet_count(); ++facet) {
            if (domain.facet_boundary_condition(facet) == 'P') {
                check(domain.has_periodic_pair(facet), "every periodic facet has a pair");
            }
        }

        SimulationConfig finfet = load_simulation_config(
            (fs::path(PHONOMC_SOURCE_DIR) / "example" / "input_finfet_stl_heat1e20.toml").string());
        finfet.grid = {2, 2, 2};
        finfet.output_folder = (root / "finfet").string();
        const SimulationDomain finfet_domain(finfet);
        check(!finfet_domain.reservoir_facets().empty(), "FinFET strict boundary configuration remains valid");
        check(!finfet_domain.is_box_geometry(), "STL geometry is distinguished from a box");
        check(finfet_domain.fast_grid_index_enabled(), "STL enables constant-time voxel lookup");
        for (const auto& center : finfet_domain.grid_centers()) {
            const int index = finfet_domain.fast_grid_index(center);
            check(index >= 0 && index < finfet_domain.grid_count(),
                  "each retained STL grid center maps to a valid grid index");
        }
        check(finfet_domain.fast_grid_index(finfet_domain.bounds_min()) >= 0,
              "an STL bounding-box cell maps to a retained active grid");

        // One deliberately reversed triangle must be repaired without changing
        // the analytical tetrahedron volume or volume-sampling correctness.
        SurfaceMesh tetra(
            {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}},
            {{0, 1, 2}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}});
        check(std::abs(tetra.volume() - (1.0 / 6.0)) <= 1e-12,
              "inconsistent tetrahedron winding is repaired before signed-volume evaluation");
        for (const auto& point : tetra.sample_volume_points(32)) {
            check(tetra.contains_point_ray_cast(point),
                  "tetrahedral volume sampler only returns points inside the closed surface");
        }
        std::mt19937_64 seeded_a(987654321);
        std::mt19937_64 seeded_b(987654321);
        const auto points_a = tetra.sample_volume_points(16, seeded_a);
        const auto points_b = tetra.sample_volume_points(16, seeded_b);
        check(points_a == points_b, "equal seeds reproduce identical volume samples");

        SimulationConfig unmatched = valid;
        unmatched.output_folder = (root / "unmatched").string();
        unmatched.boundary_position = {"relative_box", "2", "2", "2", "3", "3", "3"};
        unmatched.boundary_conditions = {BoundaryCondition::ThermalReservoir};
        unmatched.boundary_values = {300.0};
        unmatched.periodic_pair.clear();
        expect_throw([&] { (void) SimulationDomain(unmatched); }, "unmatched boundary selector is rejected");

        SimulationConfig unpaired = valid;
        unpaired.output_folder = (root / "unpaired").string();
        unpaired.boundary_position = {"relative_box", "-0.01", "0", "0", "0.01", "1", "1"};
        unpaired.boundary_conditions = {BoundaryCondition::Periodic};
        unpaired.boundary_values = {0.0};
        unpaired.periodic_pair.clear();
        expect_throw([&] { (void) SimulationDomain(unpaired); }, "unpaired periodic facet is rejected");
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
    std::cout << "All SimulationDomain tests passed\n";
    return 0;
}
