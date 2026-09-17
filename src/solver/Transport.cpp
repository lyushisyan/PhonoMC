#include "MonteCarloSolver.h"

#include "PhononMaterial.h"
#include "SimulationDomain.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

#ifdef PHONOMC_USE_OPENMP
#include <omp.h>
#endif

// 函数说明：将粒子位置映射到最近控制体网格索引。
int MonteCarloSolver::nearest_grid_index(const SimulationDomain& geometry, const Vec3& p) const {
    if (geometry.fast_grid_index_enabled()) {
        const int idx = geometry.fast_grid_index(p);
        if (idx >= 0) {
            return idx;
        }
    }
    const auto& centers = geometry.grid_centers();
    if (centers.empty()) {
        return 0;
    }
    int best = 0;
    double best_d2 = std::numeric_limits<double>::max();
    for (int i = 0; i < static_cast<int>(centers.size()); ++i) {
        const auto d = sub(p, centers[i]);
        const double d2 = dot(d, d);
        if (d2 < best_d2) {
            best_d2 = d2;
            best = i;
        }
    }
    return best;
}

// 函数说明：批量追踪粒子到下一次边界碰撞点并缓存结果。
void MonteCarloSolver::update_collision_cache(const SimulationDomain& geometry, const std::vector<int>& indices) {
    const int nidx = static_cast<int>(indices.size());
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int k = 0; k < nidx; ++k) {
        update_collision_cache_single(geometry, indices[static_cast<size_t>(k)]);
    }
    throw_if_collision_cache_failed("updating collision cache");
}

// 函数说明：为单个粒子更新下一次边界碰撞缓存。
void MonteCarloSolver::update_collision_cache_single(const SimulationDomain& geometry, int i) {
    const auto& mesh = geometry.mesh();
    if (i < 0 || i >= particles_.size()) {
        return;
    }
    const double particle_speed = norm(particles_.velocities[i]);
    if (std::isfinite(particle_speed) && particle_speed <= 1e-18) {
        particles_.collision_positions[i] = particles_.positions[i];
        particles_.collision_facets[i] = -1;
        particles_.collision_conditions[i] = 'R';
        particles_.collision_times[i] = std::numeric_limits<double>::infinity();
        particles_.collision_failed[static_cast<size_t>(i)] = static_cast<std::uint8_t>(0);
        return;
    }
    auto [cp, t_hit, fct] = geometry.trace_boundary_intersection(particles_.positions[i], particles_.velocities[i]);
    if (fct < 0 || !std::isfinite(t_hit)) {
        // Numerical correction only: try nearby interior points while keeping
        // the particle velocity unchanged. Reversing velocity here would turn
        // a geometry failure into an unphysical scattering event.
        const Vec3 original = particles_.positions[i];
        const Vec3 ext = sub(geometry.bounds_max(), geometry.bounds_min());
        const double scale = std::max({1.0, std::abs(ext[0]), std::abs(ext[1]), std::abs(ext[2])});
        const double eps = std::max(push_eps_, 1e-10 * scale);
        std::array<Vec3, 4> candidates {original, original, original, original};
        const int nearest = mesh.nearest_facet(original);
        if (nearest >= 0 && nearest < static_cast<int>(mesh.facet_normals().size())) {
            const Vec3 normal = mesh.facet_normals()[static_cast<size_t>(nearest)];
            candidates[0] = add(original, mul(normal, eps));
            candidates[1] = add(original, mul(normal, -eps));
        }
        const double speed = norm(particles_.velocities[i]);
        if (speed > 0.0 && std::isfinite(speed)) {
            const Vec3 direction = mul(particles_.velocities[i], 1.0 / speed);
            candidates[2] = add(original, mul(direction, eps));
            candidates[3] = add(original, mul(direction, -eps));
        }
        for (const Vec3& candidate : candidates) {
            if (!geometry.contains_point(candidate)) {
                continue;
            }
            std::tie(cp, t_hit, fct) = geometry.trace_boundary_intersection(candidate, particles_.velocities[i]);
            if (fct >= 0 && std::isfinite(t_hit)) {
                particles_.positions[i] = candidate;
                collision_cache_corrections_total_.fetch_add(1, std::memory_order_relaxed);
                break;
            }
        }
    }
    if (fct < 0 || !std::isfinite(t_hit)) {
        particles_.collision_positions[i] = particles_.positions[i];
        particles_.collision_facets[i] = -1;
        particles_.collision_conditions[i] = 'R';
        particles_.collision_times[i] = std::numeric_limits<double>::infinity();
        if (particles_.collision_failed[static_cast<size_t>(i)] == 0) {
            particles_.collision_failed[static_cast<size_t>(i)] = static_cast<std::uint8_t>(1);
            collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
        }
        return;
    }
    particles_.collision_failed[static_cast<size_t>(i)] = static_cast<std::uint8_t>(0);
    particles_.collision_positions[i] = cp;
    particles_.collision_facets[i] = fct;
    particles_.collision_conditions[i] = geometry.facet_boundary_condition(fct);
    particles_.collision_times[i] = std::isinf(t_hit) ? std::numeric_limits<double>::infinity() : (t_hit / time_step_);
}

void MonteCarloSolver::throw_if_collision_cache_failed(const char* context) const {
    for (size_t i = 0; i < particles_.collision_failed.size(); ++i) {
        if (particles_.collision_failed[i] == 0) {
            continue;
        }
        const Vec3& p = particles_.positions[i];
        if (particles_.collision_failed[i] == 2)
            throw std::runtime_error(
                std::string("Full-matrix boundary mode mapping failed while ") + context +
                " for particle " + std::to_string(i) + ". An inward discrete material mode and finite signed energy weight are required; no velocity-only reflection was applied.");
        throw std::runtime_error(
            std::string("Boundary intersection failed while ") + context +
            " for particle " + std::to_string(i) + " at (" +
            std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " +
            std::to_string(p[2]) + "). Check that the STL is closed and manifold.");
    }
}

void MonteCarloSolver::throw_if_excessive_collisions() const {
    const int i = excessive_collision_particle_.load(std::memory_order_relaxed);
    if (i < 0 || i >= particles_.size()) {
        return;
    }
    const Vec3& p = particles_.positions[static_cast<size_t>(i)];
    throw std::runtime_error(
        "Particle " + std::to_string(i) +
        " exceeded 64 transport events in one timestep at (" +
        std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " +
        std::to_string(p[2]) +
        "), cached facet=" + std::to_string(particles_.collision_facets[i]) +
        ". Reduce simulation.time_step or inspect narrow STL features; no recovery was applied.");
}

// 函数说明：处理边界事件（透射、吸收、周期、粗糙/镜面反射）并更新粒子态。
void MonteCarloSolver::process_boundary_collision(const SimulationDomain& geometry, int i) {
    const int facet = particles_.collision_facets[i];
    const char cond = particles_.collision_conditions[i];

    particles_.positions[i] = particles_.collision_positions[i];
    if (facet < 0 || facet >= static_cast<int>(geometry.mesh().facet_normals().size())) {
        if (uses_conservative_transport(args_)) {
            particles_.collision_failed[i] = 2;
            collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
        }
        return;
    }

    const Vec3 n = geometry.mesh().facet_normals()[facet];
    const double vn = dot(particles_.velocities[i], n);

    if (cond == 'T') {
        double absorbed_energy_ev = 0.0;
        if (phonon_ != nullptr) {
            const PhononMaterial& phonon = particle_material(i);
            constexpr double kHbarEvPs = 6.582119569e-4;
            const double neq = phonon.bose_occupation(
                background_temperature_reference(), particles_.modes[static_cast<size_t>(i)]);
            const double raw_deviation_ev = kHbarEvPs *
                phonon.mode_angular_frequency(particles_.modes[static_cast<size_t>(i)]) *
                (particles_.occupation[static_cast<size_t>(i)] - neq);
            absorbed_energy_ev = raw_deviation_ev * material_energy_weight(
                particles_.material_ids[static_cast<size_t>(i)]);
        }
#ifdef PHONOMC_USE_OPENMP
        const int worker = omp_get_thread_num();
#else
        const int worker = 0;
#endif
        boundary_ledger_.record_absorption(worker, reservoir_sampler_.reservoir_index(facet), absorbed_energy_ev);
        if (i >= 0 && i < static_cast<int>(particles_.alive.size())) {
            particles_.alive[static_cast<size_t>(i)] = static_cast<std::uint8_t>(0);
        }
        return;
    } else if (cond == 'P' && geometry.has_periodic_pair(facet)) {
        // Teleport through periodic interface with facet-pair translation.
        const Vec3 shift = geometry.periodic_shift_for_facet(facet);
        particles_.positions[i] = add(particles_.positions[i], shift);
    } else if (uses_conservative_transport(args_)) {
        // A full matrix is indexed by physical q/branch modes. Every reflection
        // must select an inward mode with its actual group velocity; changing a
        // velocity while retaining its original mode changes the transport BTE.
        const PhononMaterial& phonon = particle_material(i);
        const int rough_idx = rough_sampler_.table_index(particles_.material_ids[i], facet);
        if (rough_idx < 0) {
            particles_.collision_failed[i] = 2;
            collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
            return;
        }
        double out_occupation = particles_.occupation[i];
        const auto out_mode = select_reflected_mode(phonon, rough_idx, particles_.modes[i],
            out_occupation, particles_.occupation[i]);
        const Vec3 out_velocity = phonon.mode_group_velocity(out_mode);
        const double out_omega = phonon.mode_angular_frequency(out_mode);
        if (!std::isfinite(out_occupation) || !std::isfinite(out_omega) || !(out_omega > 0.0) ||
            !std::isfinite(norm(out_velocity)) || !(dot(out_velocity, n) < 0.0)) {
            particles_.collision_failed[i] = 2;
            collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
            return;
        }
        particles_.modes[i] = out_mode;
        particles_.velocities[i] = out_velocity;
        particles_.omega[i] = out_omega;
        particles_.occupation[i] = out_occupation;
        particles_.energies[i] = 0.0;
    } else if (cond == 'R' && phonon_ != nullptr) {
        const int material_id = particles_.material_ids[static_cast<size_t>(i)];
        const PhononMaterial& phonon = material(material_id);
        const int rough_idx = rough_sampler_.table_index(material_id, facet);
        if (rough_idx >= 0) {
            const double in_occ = particles_.occupation[i];
            double out_occ = in_occ;
            const std::array<int, 2> out_mode =
                select_reflected_mode(phonon, rough_idx, particles_.modes[i], out_occ, in_occ);
            particles_.modes[i] = out_mode;
            particles_.velocities[i] = phonon.mode_group_velocity(particles_.modes[i]);
            if (dot(particles_.velocities[i], n) > 0.0) {
                particles_.velocities[i] = sub(particles_.velocities[i], mul(n, 2.0 * dot(particles_.velocities[i], n)));
            }
            particles_.omega[i] = phonon.mode_angular_frequency(particles_.modes[i]);
            particles_.occupation[i] = out_occ;
            particles_.energies[i] = 0.0;
        } else {
            rough_events_total_.fetch_add(1, std::memory_order_relaxed);
            rough_fallback_missing_rough_data_.fetch_add(1, std::memory_order_relaxed);
            const double p_spec = compute_roughness_specularity(geometry, phonon, i, facet);
            std::uniform_real_distribution<double> U(0.0, 1.0);
            auto& rng = thread_rng();
            if (U(rng) <= p_spec) {
                rough_specular_selected_.fetch_add(1, std::memory_order_relaxed);
                particles_.velocities[i] = sub(particles_.velocities[i], mul(n, 2.0 * vn));
            } else {
                rough_diffuse_selected_.fetch_add(1, std::memory_order_relaxed);
                rough_fallback_global_random_.fetch_add(1, std::memory_order_relaxed);
                // Missing rough-mode data must not thermalize an adiabatic
                // wall. Keep mode and occupation, randomizing direction only.
                const double speed = std::max(1e-9, norm(particles_.velocities[i]));
                Vec3 dir = random_unit_vector();
                if (dot(dir, n) > 0.0) {
                    dir = mul(dir, -1.0);
                }
                particles_.velocities[i] = mul(dir, speed);
                particles_.omega[i] = phonon.mode_angular_frequency(particles_.modes[i]);
            }
        }
    } else {
        // Unmatched periodic and generic boundaries: specular reflection.
        particles_.velocities[i] = sub(particles_.velocities[i], mul(n, 2.0 * vn));
    }

    particles_.positions[i] = add(particles_.positions[i], mul(particles_.velocities[i], push_eps_ / std::max(norm(particles_.velocities[i]), 1e-12)));
    particles_.grid_ids[i] = nearest_grid_index(geometry, particles_.positions[i]);
    const int destination_material = material_index_at(particles_.positions[i]);
    if (uses_conservative_transport(args_) && destination_material != particles_.material_ids[i]) {
        particles_.collision_failed[i] = 2;
        collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    particles_.material_ids[i] = destination_material;
}

// 函数说明：将数值误差导致越界的粒子重采样回几何体内部，并重建碰撞缓存。
void MonteCarloSolver::throw_if_particles_escaped(const SimulationDomain& geometry) const {
    const auto& bmin = geometry.bounds_min();
    const auto& bmax = geometry.bounds_max();
    const Vec3 ext = sub(bmax, bmin);
    const double scale = std::max({1.0, std::abs(ext[0]), std::abs(ext[1]), std::abs(ext[2])});
    const double tol = 1e-10 * scale;
    std::atomic<int> escaped {-1};
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particles_.size(); ++i) {
        const Vec3& p = particles_.positions[static_cast<size_t>(i)];
        bool outside = false;
        for (int axis = 0; axis < 3; ++axis)
            outside = outside || !std::isfinite(p[axis]) ||
                p[axis] < bmin[axis] - tol || p[axis] > bmax[axis] + tol;
        if (!outside && !geometry.is_box_geometry()) outside = !geometry.contains_point(p);
        if (outside) {
            int expected = -1;
            escaped.compare_exchange_strong(expected, i, std::memory_order_relaxed);
        }
    }
    const int i = escaped.load(std::memory_order_relaxed);
    if (i < 0) return;
    const Vec3& p = particles_.positions[static_cast<size_t>(i)];
    throw std::runtime_error("Particle " + std::to_string(i) + " escaped the domain at (" +
        std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) +
        "), cached facet=" + std::to_string(particles_.collision_facets[i]) +
        ". Refusing to relocate or rethermalize it; inspect the geometry and time step.");
}

// 函数说明：在单个时间步内推进粒子运动并处理可能的多次边界碰撞。
void MonteCarloSolver::advance_particle(const SimulationDomain& geometry, int i, double dt_remaining) {
    if (i < 0 || i >= particles_.size()) {
        return;
    }
    if (i < static_cast<int>(particles_.alive.size()) && particles_.alive[static_cast<size_t>(i)] == 0) {
        return;
    }
    double remaining = dt_remaining;
    int guard = 0;
    while (remaining > 1e-14 && guard < 64) {
        ++guard;
        const double t_hit = std::isinf(particles_.collision_times[i])
            ? std::numeric_limits<double>::infinity()
            : particles_.collision_times[i] * time_step_;
        int next_material = -1;
        const double t_interface = next_interface_time(i, next_material);
        const double t_event = std::min(t_hit, t_interface);
        if (!std::isfinite(t_event) || t_event > remaining) {
            particles_.positions[i] = add(particles_.positions[i], mul(particles_.velocities[i], remaining));
            if (std::isfinite(particles_.collision_times[i])) {
                particles_.collision_times[i] -= remaining / time_step_;
            }
            remaining = 0.0;
        } else if (t_interface < t_hit) {
            particles_.positions[i] = add(
                particles_.positions[i], mul(particles_.velocities[i], t_interface));
            remaining -= std::max(0.0, t_interface);
            process_material_interface(i, next_material);
            if (particles_.collision_failed[static_cast<size_t>(i)] != 0) break;
            const double speed = std::max(norm(particles_.velocities[i]), 1e-12);
            particles_.positions[i] = add(
                particles_.positions[i], mul(particles_.velocities[i], push_eps_ / speed));
            particles_.grid_ids[i] = nearest_grid_index(geometry, particles_.positions[i]);
            update_collision_cache_single(geometry, i);
            if (particles_.collision_failed[static_cast<size_t>(i)] != 0) {
                break;
            }
        } else {
            // Move to collision and process boundary event.
            particles_.positions[i] = particles_.collision_positions[i];
            remaining -= std::max(0.0, t_hit);
            process_boundary_collision(geometry, i);
            if (particles_.alive[static_cast<size_t>(i)] == 0 || particles_.collision_failed[static_cast<size_t>(i)] != 0) {
                break;
            }
            update_collision_cache_single(geometry, i);
            if (particles_.collision_failed[static_cast<size_t>(i)] != 0) {
                break;
            }
            if (particles_.collision_times[i] * time_step_ < 1e-12) {
                particles_.positions[i] = add(particles_.positions[i], mul(particles_.velocities[i], 1e-12));
                update_collision_cache_single(geometry, i);
                if (particles_.collision_failed[static_cast<size_t>(i)] != 0) {
                    break;
                }
            }
        }
    }
    if (remaining > 1e-14 && particles_.alive[static_cast<size_t>(i)] != 0 &&
        particles_.collision_failed[static_cast<size_t>(i)] == 0) {
        // Report outside the OpenMP region; do not teleport the particle or
        // discard its remaining flight time and continue a corrupted trajectory.
        int expected = -1;
        excessive_collision_particle_.compare_exchange_strong(expected, i, std::memory_order_relaxed);
    }
}
