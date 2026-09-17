#include "MonteCarloSolver.h"

#include "PhononMaterial.h"
#include "SimulationDomain.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

void MonteCarloSolver::initialize_material_layout(const SimulationDomain& geometry) {
    grid_.material_ids.resize(geometry.grid_centers().size(), 0);
    for (size_t grid = 0; grid < geometry.grid_centers().size(); ++grid) {
        grid_.material_ids[grid] = material_index_at(geometry.grid_centers()[grid]);
    }
    if (materials_.size() <= 1) {
        return;
    }
    const int cells = (args_.material_interface_axis == 0) ? args_.grid.nx
        : (args_.material_interface_axis == 1) ? args_.grid.ny : args_.grid.nz;
    const double lo = geometry.bounds_min()[static_cast<size_t>(args_.material_interface_axis)];
    const double hi = geometry.bounds_max()[static_cast<size_t>(args_.material_interface_axis)];
    const double cell_width = (hi - lo) / static_cast<double>(cells);
    const double tolerance = std::max(1e-10, 1e-10 * std::abs(hi - lo));
    for (double plane : args_.material_interface_positions) {
        if (!std::isfinite(plane) || !(plane > lo) || !(plane < hi)) {
            throw std::runtime_error("Each material interface must lie strictly inside the loaded geometry bounds.");
        }
        const double cell_coordinate = (plane - lo) / cell_width;
        if (std::abs(cell_coordinate - std::round(cell_coordinate)) > tolerance / cell_width) {
            throw std::runtime_error(
                "Each materials.interface_positions value must coincide with a grid-cell boundary.");
        }
    }
    std::vector<bool> occupied(materials_.size(), false);
    for (int id : grid_.material_ids) occupied.at(static_cast<size_t>(id)) = true;
    if (std::find(occupied.begin(), occupied.end(), false) != occupied.end()) {
        throw std::runtime_error("Every material layer must contain at least one active geometry cell.");
    }
    const auto& normals = geometry.mesh().facet_normals();
    const int axis = args_.material_interface_axis;
    for (int facet = 0; facet < geometry.mesh().facet_count(); ++facet) {
        const char condition = geometry.facet_boundary_condition(facet);
        if (condition == 'P' &&
            std::abs(normals[static_cast<size_t>(facet)][static_cast<size_t>(axis)]) > 0.99) {
            throw std::runtime_error(
                "Periodic boundaries normal to materials.interface_axis are not supported, "
                "because they would join different end materials.");
        }
        if (condition == 'T' &&
            std::abs(normals[static_cast<size_t>(facet)][static_cast<size_t>(axis)]) < 0.99) {
            throw std::runtime_error(
                "In the layered-material model, thermal-reservoir facets must be normal to "
                "materials.interface_axis so each reservoir touches one material.");
        }
    }
    std::cout << "[init] layered materials=" << materials_.size()
              << ", interface_axis=" << "xyz"[args_.material_interface_axis]
              << ", interface_model=" << to_string(args_.material_interface_model) << '\n';
}

// 函数说明：初始化粒子主状态与碰撞缓存，建立时间推进的初始条件。
void MonteCarloSolver::initialize_particles(const SimulationDomain& geometry) {
    using Clock = std::chrono::steady_clock;
    auto step_begin = Clock::now();
    auto begin_step = [&](int idx, int total, const std::string& name) {
        step_begin = Clock::now();
        std::cout << "[init] " << idx << "/" << total << " " << name << "...\n";
    };
    auto end_step = [&]() {
        const auto step_end = Clock::now();
        const double sec = std::chrono::duration<double>(step_end - step_begin).count();
        std::cout << "        done (" << std::fixed << std::setprecision(2) << sec << " s)\n";
    };
    constexpr int total_steps = 8;

    const auto& mesh = geometry.mesh();
    begin_step(1, total_steps, "Sampling particle positions and assigning initial grid IDs");
    auto positions = geometry.sample_volume_points(particles_.size(), rng_);
    if (materials_.size() > 1) {
        const double smallest_volume = *std::min_element(particles_.material_spatial_volumes_a3.begin(),
                                                        particles_.material_spatial_volumes_a3.end());
        std::uniform_real_distribution<double> uniform(0.0,1.0);
        // Rejection from the actual geometry also handles clipped STL cells.
        size_t accepted=0;
        while (accepted < positions.size()) {
            const auto candidates=geometry.sample_volume_points(static_cast<int>(positions.size()-accepted),rng_);
            for (const auto& position:candidates) {
                const int mid=material_index_at(position);
                if (uniform(rng_) < smallest_volume/particles_.spatial_volume(mid,particle_spatial_weight_a3_))
                    positions[accepted++]=position;
            }
        }
    }
    std::copy(positions.begin(), positions.end(), particles_.positions.begin());
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particles_.size(); ++i) {
        particles_.grid_ids[i] = nearest_grid_index(geometry, particles_.positions[i]);
        particles_.material_ids[i] = material_index_at(particles_.positions[i]);
    }
    end_step();

    begin_step(2, total_steps, "Assigning phonon modes");
    initialize_particle_modes();
    end_step();
    begin_step(3, total_steps, "Initializing particle temperatures");
    initialize_particle_temperatures(geometry);
    end_step();
    begin_step(4, total_steps, "Initializing particle velocities");
    initialize_particle_velocities();
    end_step();
    begin_step(5, total_steps, "Building reservoir injection tables");
    initialize_reservoir_injection(geometry);
    end_step();
    begin_step(6, total_steps, "Precomputing rough-boundary scattering tables");
    initialize_rough_boundary_scattering(geometry);
    end_step();
    begin_step(7, total_steps, "Initializing particle state arrays and collision cache");
    const double tmin = particles_.temperatures.empty() ? 300.0 : *std::min_element(particles_.temperatures.begin(), particles_.temperatures.end());
    const int nsv = std::max(1, geometry.grid_count());
    grid_.temperatures.assign(static_cast<size_t>(nsv), tmin);
    grid_.energy_density.assign(static_cast<size_t>(nsv), 0.0);
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particles_.size(); ++i) {
        const PhononMaterial& phonon = particle_material(i);
        particles_.omega[i] = phonon.mode_angular_frequency(particles_.modes[i]);
        particles_.occupation[i] = uses_linearized_transport(args_)
            ? phonon.linearized_bose_occupation(particles_.temperatures[i], background_temperature_, particles_.modes[i])
            : phonon.bose_occupation(particles_.temperatures[i], particles_.modes[i]);
        particles_.energies[i] = 0.0;
    }
    apply_initial_mode_excitation();

    std::vector<int> all_idx(static_cast<size_t>(particles_.size()));
    std::iota(all_idx.begin(), all_idx.end(), 0);
    update_collision_cache(geometry, all_idx);
    end_step();
    begin_step(8, total_steps, "Computing initial temperature field from particle energy");
    update_particle_temperatures(geometry);
    end_step();
    assert(particles_.aligned());
}

// 函数说明：按材料活跃模态集合为每个粒子分配初始声子模态。
void MonteCarloSolver::initialize_particle_modes() {
    for (int i = 0; i < particles_.size(); ++i) {
        particles_.modes[i] = particle_material(i).sample_active_mode(rng_);
    }
}

void MonteCarloSolver::apply_initial_mode_excitation() {
    if (!args_.mode_excitation_enabled ||
        !(args_.mode_excitation_occupation_multiplier > 0.0) ||
        std::abs(args_.mode_excitation_occupation_multiplier - 1.0) < 1e-15) {
        return;
    }
    constexpr double kTwoPi = 6.283185307179586476925286766559;
    long long modified = 0;
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for reduction(+:modified)
#endif
    for (int i = 0; i < particles_.size(); ++i) {
        const double frequency_thz = particles_.omega[static_cast<size_t>(i)] / kTwoPi;
        if (!std::isfinite(frequency_thz)) {
            continue;
        }
        if (args_.mode_excitation_frequency_min >= 0.0 &&
            frequency_thz < args_.mode_excitation_frequency_min) {
            continue;
        }
        if (args_.mode_excitation_frequency_max >= 0.0 &&
            frequency_thz > args_.mode_excitation_frequency_max) {
            continue;
        }
        particles_.occupation[static_cast<size_t>(i)] *=
            args_.mode_excitation_occupation_multiplier;
        ++modified;
    }
    std::cout << "[init] mode_excitation enabled: frequency_min="
              << args_.mode_excitation_frequency_min << " THz"
              << ", frequency_max=" << args_.mode_excitation_frequency_max << " THz"
              << ", occupation_multiplier=" << args_.mode_excitation_occupation_multiplier
              << ", modified_particles=" << modified << '\n';
}

// 函数说明：依据边界温度与初始策略设置粒子温度场。
void MonteCarloSolver::initialize_particle_temperatures(const SimulationDomain& geometry) {
    std::fill(particles_.temperatures.begin(), particles_.temperatures.end(), 300.0);
    constexpr double kNoReservoirCold = 299.0;
    constexpr double kNoReservoirHot = 301.0;
    std::vector<double> finite_res_vals;
    finite_res_vals.reserve(geometry.reservoir_values().size());
    for (double v : geometry.reservoir_values()) {
        if (std::isfinite(v)) {
            finite_res_vals.push_back(v);
        }
    }

    double tmin = kNoReservoirCold;
    double tmax = kNoReservoirHot;
    if (!finite_res_vals.empty()) {
        tmin = *std::min_element(finite_res_vals.begin(), finite_res_vals.end());
        tmax = *std::max_element(finite_res_vals.begin(), finite_res_vals.end());
    }
    if (tmin > tmax) {
        std::swap(tmin, tmax);
    }
    const double tmean = 0.5 * (tmin + tmax);
    if (args_.initial_temperature.mode == InitialTemperatureMode::Uniform) {
        const double uniform_temperature = args_.initial_temperature.uniform_temperature;
        std::fill(particles_.temperatures.begin(), particles_.temperatures.end(), uniform_temperature);
        std::cout << "[init] initial_temperature mode=uniform, T0=" << uniform_temperature << " K\n";
    } else {
        const auto& rf = geometry.reservoir_facets();
        const auto& centroids = geometry.mesh().facet_centroids();
        int cold_facet = -1;
        int hot_facet = -1;
        double cold_t = std::numeric_limits<double>::infinity();
        double hot_t = -std::numeric_limits<double>::infinity();
        for (int facet : rf) {
            if (facet < 0 || facet >= static_cast<int>(centroids.size())) {
                continue;
            }
            const double Tres = geometry.reservoir_value_for_facet(facet, std::numeric_limits<double>::quiet_NaN());
            if (!std::isfinite(Tres)) {
                continue;
            }
            if (Tres < cold_t) {
                cold_t = Tres;
                cold_facet = facet;
            }
            if (Tres > hot_t) {
                hot_t = Tres;
                hot_facet = facet;
            }
        }

        const bool valid_pair = (cold_facet >= 0 && hot_facet >= 0 && hot_t > cold_t + 1e-12);
        if (!valid_pair) {
            std::fill(particles_.temperatures.begin(), particles_.temperatures.end(), tmean);
            std::cout << "[init] initial_temperature mode=linear, fallback_to_uniform_mean="
                      << tmean << " K (insufficient hot/cold reservoirs)\n";
            return;
        }

        const Vec3 p_cold = centroids[static_cast<size_t>(cold_facet)];
        const Vec3 p_hot = centroids[static_cast<size_t>(hot_facet)];
        const Vec3 axis = sub(p_hot, p_cold);
        const double axis_len2 = dot(axis, axis);
        if (!(axis_len2 > 1e-20)) {
            std::fill(particles_.temperatures.begin(), particles_.temperatures.end(), tmean);
            std::cout << "[init] initial_temperature mode=linear, fallback_to_uniform_mean="
                      << tmean << " K (degenerate hot/cold direction)\n";
            return;
        }

#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < particles_.size(); ++i) {
            const Vec3 rel = sub(particles_.positions[i], p_cold);
            const double s = std::clamp(dot(rel, axis) / axis_len2, 0.0, 1.0);
            particles_.temperatures[static_cast<size_t>(i)] = cold_t + s * (hot_t - cold_t);
        }
        std::cout << "[init] initial_temperature mode=linear, cold=" << cold_t
                  << " K (facet " << cold_facet << "), hot=" << hot_t
                  << " K (facet " << hot_facet << ")\n";
    }
}

// 函数说明：采样随机单位方向，用于漫反射与方向随机化过程。
MonteCarloSolver::Vec3 MonteCarloSolver::random_unit_vector() {
    auto& rng = thread_rng();
    std::normal_distribution<double> N(0.0, 1.0);
    Vec3 v {N(rng), N(rng), N(rng)};
    const double n = norm(v);
    if (n <= 1e-12) {
        return {1.0, 0.0, 0.0};
    }
    return mul(v, 1.0 / n);
}

// 函数说明：根据粒子模态查询群速度并写入速度场。
void MonteCarloSolver::initialize_particle_velocities() {
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particles_.size(); ++i) {
        particles_.velocities[i] = particle_material(i).mode_group_velocity(particles_.modes[i]);
    }
}
