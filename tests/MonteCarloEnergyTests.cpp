#include "MonteCarloSolver.h"
#include "FluxWeightedWindow.h"
#include "PhononMaterial.h"
#include "SimulationConfig.h"
#include "SimulationDomain.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
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

void check_close(double actual, double expected, double relative_tolerance, const std::string& message) {
    const double scale = std::max({1e-14, std::abs(actual), std::abs(expected)});
    check(std::abs(actual - expected) <= relative_tolerance * scale,
          message + " (actual=" + std::to_string(actual) +
          ", expected=" + std::to_string(expected) + ")");
}
}  // namespace

int main() {
    {
        const std::vector<double> outgoing_flux {4.0, 6.0};
        const std::vector<double> incoming_flux {5.0, 5.0};
        const std::vector<double> specularity {0.5, 0.2};
        const std::vector<int> specular_match {0, 1};
        const std::vector<double> balanced_residual =
            phonomc_detail::residual_diffuse_creation_rates(
                outgoing_flux, incoming_flux, specularity, specular_match);
        check_close(balanced_residual[0], 1.5, 1e-14,
                    "diffuse residual removes the first specular channel");
        check_close(balanced_residual[1], 5.0, 1e-14,
                    "diffuse residual removes the second specular channel");
        const double diffuse_incoming =
            incoming_flux[0] * (1.0 - specularity[0]) +
            incoming_flux[1] * (1.0 - specularity[1]);
        check_close(balanced_residual[0] + balanced_residual[1], diffuse_incoming, 1e-14,
                    "specular plus diffuse channels satisfy shell-wise detailed balance");

        const std::vector<double> residual_rates {0.5, -1.0, 2.0, 1.5};
        const std::vector<int> frequency_order {2, 0, 3, 1};
        const std::vector<double> residual_prefix =
            phonomc_detail::ordered_nonnegative_prefix(residual_rates, frequency_order);
        check(residual_prefix.size() == 5,
              "residual-flux prefix has one sentinel entry");
        check_close(residual_prefix[1], 2.0, 1e-14,
                    "residual-flux prefix follows frequency-sorted mode order");
        check_close(residual_prefix[3], 4.0, 1e-14,
                    "residual-flux prefix accumulates nonnegative creation rates");
        check_close(residual_prefix[4], 4.0, 1e-14,
                    "residual-flux prefix rejects negative creation rates");
        check(phonomc_detail::flux_weighted_window_index(residual_prefix, 0, 3, 0.49) == 0,
              "residual sampler uses the first mode below its cumulative boundary");
        check(phonomc_detail::flux_weighted_window_index(residual_prefix, 0, 3, 0.51) == 1,
              "residual sampler selects from the specular-depleted outgoing flux");

        const std::vector<double> prefix {0.0, 0.01, 1.0};
        check(phonomc_detail::flux_weighted_window_index(prefix, 0, 2, 0.005) == 0,
              "flux sampler selects the one-percent interval below its probability boundary");
        check(phonomc_detail::flux_weighted_window_index(prefix, 0, 2, 0.02) == 1,
              "flux sampler assigns the remaining probability to the high-flux interval");
        const std::vector<double> leading_zero {0.0, 0.0, 1.0};
        check(phonomc_detail::flux_weighted_window_index(leading_zero, 0, 2, 0.0) == 1,
              "flux sampler never selects a zero-normal-flux interval");
        const std::vector<double> subwindow {0.0, 0.1, 0.3, 1.0};
        check(phonomc_detail::flux_weighted_window_index(subwindow, 1, 3, 0.0) == 1,
              "flux sampler respects an elastic subwindow begin index");
        check(phonomc_detail::flux_weighted_window_index(subwindow, 1, 3, 0.3) == 2,
              "flux sampler respects weights normalized inside an elastic subwindow");
        const std::vector<double> zero_window {0.0, 0.0, 0.0};
        check(phonomc_detail::flux_weighted_window_index(zero_window, 0, 2, 0.5) == -1,
              "flux sampler rejects a zero-normal-flux window");
    }

    const auto tag = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const fs::path root = fs::temp_directory_path() / ("phonomc-energy-tests-" + std::to_string(tag));
    fs::create_directories(root);
    const fs::path input = root / "energy.toml";
    {
        std::ofstream out(input);
        out << R"([geometry]
model = "box"
sizes = [10, 10, 10]

[simulation]
particle_count = 4096
time_step = 0.05
iterations = 20
convergence_write_interval = 20
random_seed = 24680
compute_kappa = false
initial_temperature = 300
background_temperature = 300
lifetime_temperature = "local"
grid_xyz = [2, 2, 2]

[boundary]
boundary_position = [
  [-0.01, 0, 0, 0.01, 1, 1],
  [0.99, 0, 0, 1.01, 1, 1],
  [0, -0.01, 0, 1, 0.01, 1],
  [0, 0.99, 0, 1, 1.01, 1],
  [0, 0, -0.01, 1, 1, 0.01],
  [0, 0, 0.99, 1, 1, 1.01]
]
boundary_conditions = ["P", "P", "P", "P", "P", "P"]
boundary_values = [0, 0, 0, 0, 0, 0]
periodic_pair = [
  [-0.01, 0, 0, 0.01, 1, 1],
  [0.99, 0, 0, 1.01, 1, 1],
  [0, -0.01, 0, 1, 0.01, 1],
  [0, 0.99, 0, 1, 1.01, 1],
  [0, 0, -0.01, 1, 1, 0.01],
  [0, 0, 0.99, 1, 1, 1.01]
]

[heat_source]
enabled = true
profile = "uniform"
min = [0, 0, 0]
max = [10, 10, 10]
power_density = 1e15
time_profile = "constant"

[io]
material_folder = "missing-material"
output_folder = ""
)";
    }

#if defined(_WIN32)
    _putenv_s("PHONOMC_ALLOW_SYNTHETIC_MATERIAL", "1");
#else
    setenv("PHONOMC_ALLOW_SYNTHETIC_MATERIAL", "1", 1);
#endif

    try {
        const SimulationConfig config = load_simulation_config(input.string());
        const SimulationDomain domain(config);
        const PhononMaterial material(config, 0);
        MonteCarloSolver solver(config, domain, material);
        const double initial_energy_ev = solver.total_thermal_energy_ev();
        for (int step = 0; step < config.iterations; ++step) {
            solver.run_timestep();
        }

        constexpr double wm3_to_evpsa3 = 6.241509074e-24;
        const double expected_source_energy_ev = config.heat_source_power_density * wm3_to_evpsa3 *
            domain.volume() * config.time_step * static_cast<double>(config.iterations);
        check(solver.particle_count() == static_cast<int>(config.particle_count),
              "occupation heat source keeps the carrier count fixed in a periodic domain");
        check_close(solver.total_heat_source_energy_ev(), expected_source_energy_ev, 1e-11,
                    "integrated occupation source equals prescribed physical energy");
        check_close(solver.total_thermal_energy_ev() - initial_energy_ev, expected_source_energy_ev, 2e-7,
                    "periodic source-only domain gains exactly the prescribed thermal energy");
        check(std::abs(solver.total_lifetime_energy_residual_ev()) <= 1e-9,
              "finite-step lifetime scattering conserves represented energy");
        check(std::abs(solver.total_energy_balance_residual_ev()) <= 1e-8,
              "fixed background makes occupation-source energy close globally");

        SimulationConfig full_phonon_config = config;
        full_phonon_config.background_temperature = 0.0;
        full_phonon_config.iterations = 5;
        const SimulationDomain full_phonon_domain(full_phonon_config);
        MonteCarloSolver full_phonon_solver(full_phonon_config, full_phonon_domain, material);
        const double full_phonon_initial_energy_ev = full_phonon_solver.total_thermal_energy_ev();
        for (int step = 0; step < full_phonon_config.iterations; ++step) {
            full_phonon_solver.run_timestep();
        }
        const double expected_full_phonon_energy_ev =
            full_phonon_config.heat_source_power_density * wm3_to_evpsa3 *
            full_phonon_domain.volume() * full_phonon_config.time_step *
            static_cast<double>(full_phonon_config.iterations);
        check_close(
            full_phonon_solver.total_thermal_energy_ev() - full_phonon_initial_energy_ev,
            expected_full_phonon_energy_ev,
            2e-7,
            "zero-K background full-phonon mode preserves source energy");
        check(std::abs(full_phonon_solver.total_energy_balance_residual_ev()) <= 1e-8,
              "zero-K background full-phonon mode closes the energy ledger");

        SimulationConfig square_config = config;
        square_config.lifetime_temperature_is_local = false;
        square_config.lifetime_temperature = 300.0;
        square_config.time_step = 0.07;
        square_config.iterations = 20;
        square_config.heat_source_time_profile = HeatSourceTimeProfile::Square;
        square_config.heat_source_period = 0.2;
        square_config.heat_source_on_duration = 0.05;
        const SimulationDomain square_domain(square_config);
        MonteCarloSolver square_solver(square_config, square_domain, material);
        const double square_initial_energy_ev = square_solver.total_thermal_energy_ev();
        for (int step = 0; step < square_config.iterations; ++step) {
            square_solver.run_timestep();
        }
        // 1.4 ps contains exactly seven periods, each with 0.05 ps on-time.
        const double expected_square_energy_ev = square_config.heat_source_power_density * wm3_to_evpsa3 *
            square_domain.volume() * (7.0 * 0.05);
        check_close(square_solver.total_heat_source_energy_ev(), expected_square_energy_ev, 1e-11,
                    "square source integrates transitions inside a timestep exactly");
        check_close(square_solver.total_thermal_energy_ev() - square_initial_energy_ev,
                    expected_square_energy_ev, 2e-7,
                    "square occupation source preserves total energy in a periodic domain");
        check(std::abs(square_solver.total_energy_balance_residual_ev()) <= 1e-8,
              "square source closes the step-wise energy ledger");

        SimulationConfig rough_config = config;
        rough_config.time_step = 0.1;
        rough_config.iterations = 50;
        rough_config.boundary_conditions.assign(6, BoundaryCondition::Rough);
        rough_config.boundary_values.assign(6, 1.0);
        rough_config.periodic_pair.clear();
        const SimulationDomain rough_domain(rough_config);
        MonteCarloSolver rough_solver(rough_config, rough_domain, material);
        const double rough_initial_energy_ev = rough_solver.total_thermal_energy_ev();
        for (int step = 0; step < rough_config.iterations; ++step) {
            rough_solver.run_timestep();
        }
        const double expected_rough_energy_ev = rough_config.heat_source_power_density * wm3_to_evpsa3 *
            rough_domain.volume() * rough_config.time_step * static_cast<double>(rough_config.iterations);
        check(rough_solver.particle_count() == static_cast<int>(rough_config.particle_count),
              "adiabatic rough boundaries preserve the carrier count");
        check_close(rough_solver.total_thermal_energy_ev() - rough_initial_energy_ev,
                    expected_rough_energy_ev, 5e-7,
                    "adiabatic rough scattering preserves occupation-source energy");
        check(std::abs(rough_solver.total_energy_balance_residual_ev()) <= 1e-8,
              "rough-boundary occupation source closes the energy ledger");
    } catch (const std::exception& ex) {
        ++failures;
        std::cerr << "Unexpected exception: " << ex.what() << '\n';
    }

    std::error_code ec;
    fs::remove_all(root, ec);
    if (failures != 0) {
        std::cerr << failures << " energy test assertion(s) failed\n";
        return 1;
    }
    std::cout << "Monte Carlo energy conservation tests passed\n";
    return 0;
}
