#include "solver/TransportOperator.h"

#include "PhononMaterial.h"
#include "SimulationDomain.h"
#include "solver/BoundaryLedger.h"
#include "solver/GridLocation.h"
#include "solver/InterfaceInteraction.h"
#include "solver/LayeredMaterialLayout.h"
#include "solver/ReservoirSampler.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#ifdef PHONOMC_USE_OPENMP
#include <omp.h>
#endif

namespace phonomc {
namespace {
Vec3 add(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 sub(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 mul(const Vec3& a, double s) {
    return {a[0] * s, a[1] * s, a[2] * s};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const Vec3& a) {
    return std::sqrt(dot(a, a));
}

int worker_index() {
#ifdef PHONOMC_USE_OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}
} // namespace

std::mt19937_64& TransportRandomStreams::for_worker(int worker) const {
#ifdef PHONOMC_USE_OPENMP
    return workers[static_cast<size_t>(worker)];
#else
    (void)worker;
    return serial;
#endif
}

void TransportRandomStreams::validate() const {
#ifdef PHONOMC_USE_OPENMP
    if (workers.size() < static_cast<size_t>(std::max(1, omp_get_max_threads()))) {
        throw std::invalid_argument("Transport requires an RNG stream for every OpenMP worker.");
    }
#endif
}

TransportOperator::TransportOperator(TransportOptions options,
                                     const SimulationDomain& geometry,
                                     const LayeredMaterialLayout& layout,
                                     const std::vector<const PhononMaterial*>& materials,
                                     std::vector<double> energy_weights,
                                     const RoughBoundarySampler& rough,
                                     const InterfaceInteraction& interfaces,
                                     const ReservoirSampler& reservoirs)
    : options_(options), geometry_(geometry), layout_(layout), materials_(materials),
      energy_weights_(std::move(energy_weights)), rough_sampler_(rough),
      interface_interaction_(interfaces), reservoir_sampler_(reservoirs) {
    if (!std::isfinite(options_.time_step_ps) || options_.time_step_ps <= 0 ||
        !std::isfinite(options_.push_distance_a) || options_.push_distance_a < 0 ||
        !std::isfinite(options_.boundary.reference_temperature) ||
        options_.boundary.reference_temperature < 0 || materials_.empty() ||
        energy_weights_.size() != materials_.size()) {
        throw std::invalid_argument("Invalid transport options or material weights.");
    }
    for (size_t i = 0; i < materials_.size(); ++i) {
        if (!materials_[i] || !std::isfinite(energy_weights_[i]) || energy_weights_[i] <= 0) {
            throw std::invalid_argument("Invalid transport material or energy weight.");
        }
    }
}

void TransportOperator::ensure_usable() const {
    if (failed_) {
        throw std::runtime_error(
            "Cannot resume a failed transport operator; construct a new simulation.");
    }
}

void TransportOperator::begin_step() {
    ensure_usable();
    excessive_collision_particle_.store(-1, std::memory_order_relaxed);
}

void TransportOperator::finish_advance(const ParticleStorage& particles, const char* context) {
    throw_if_excessive_collisions(particles);
    throw_if_collision_cache_failed(particles, context);
}

void TransportOperator::advance_existing(ParticleStorage& particles,
                                         BoundaryLedger& ledger,
                                         TransportRandomStreams rngs) {
    ensure_usable();
    rngs.validate();
    assert(particles.aligned());
    const int n_before = particles.size();
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int i = 0; i < n_before; ++i) {
        const int worker = worker_index();
        advance_particle(
            particles, i, options_.time_step_ps, ledger, worker, rngs.for_worker(worker));
        if (particles.alive[static_cast<size_t>(i)] != 0) {
            particles.grid_ids[i] = nearest_grid_index(geometry_, particles.positions[i]);
        }
    }
    finish_advance(particles, "advancing particles");
}

void TransportOperator::advance_injected(ParticleStorage& particles,
                                         const std::vector<std::pair<int, double>>& injected,
                                         BoundaryLedger& ledger,
                                         TransportRandomStreams rngs) {
    ensure_usable();
    rngs.validate();
    assert(particles.aligned());
    for (const auto& [index, emission_time] : injected) {
        if (index < 0 || index >= particles.size() || !std::isfinite(emission_time) ||
            emission_time < 0 || emission_time > options_.time_step_ps) {
            throw std::invalid_argument("Injected carriers require valid indices and emission times within the timestep.");
        }
    }
    const int ninj = static_cast<int>(injected.size());
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int k = 0; k < ninj; ++k) {
        const auto [idx, dt_in] = injected[static_cast<size_t>(k)];
        const double remain = options_.time_step_ps - dt_in;
        if (remain > 1e-14) {
            const int worker = worker_index();
            advance_particle(particles, idx, remain, ledger, worker, rngs.for_worker(worker));
        }
        if (particles.alive[static_cast<size_t>(idx)] != 0) {
            particles.grid_ids[idx] = nearest_grid_index(geometry_, particles.positions[idx]);
        }
    }
    finish_advance(particles, "advancing injected particles");
}

// 函数说明：批量追踪粒子到下一次边界碰撞点并缓存结果。
void TransportOperator::update_collision_cache(ParticleStorage& particles,
                                               const std::vector<int>& indices) {
    ensure_usable();
    assert(particles.aligned());
    for (int index : indices) {
        if (index < 0 || index >= particles.size()) {
            throw std::invalid_argument("Collision-cache batch contains an invalid particle index.");
        }
    }
    const int nidx = static_cast<int>(indices.size());
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int k = 0; k < nidx; ++k) {
        update_collision_cache_single(particles, indices[static_cast<size_t>(k)]);
    }
    throw_if_collision_cache_failed(particles, "updating collision cache");
}

// 函数说明：为单个粒子更新下一次边界碰撞缓存。
void TransportOperator::update_collision_cache_single(ParticleStorage& particles, int i) {
    const auto& geometry = geometry_;
    const auto& mesh = geometry.mesh();
    const double particle_speed = norm(particles.velocities[i]);
    if (std::isfinite(particle_speed) && particle_speed <= 1e-18) {
        particles.collision_positions[i] = particles.positions[i];
        particles.collision_facets[i] = -1;
        particles.collision_conditions[i] = 'R';
        particles.collision_times[i] = std::numeric_limits<double>::infinity();
        particles.collision_failed[static_cast<size_t>(i)] = static_cast<std::uint8_t>(0);
        return;
    }
    auto [cp, t_hit, fct] =
        geometry.trace_boundary_intersection(particles.positions[i], particles.velocities[i]);
    if (fct < 0 || !std::isfinite(t_hit)) {
        // Numerical correction only: try nearby interior points while keeping
        // the particle velocity unchanged. Reversing velocity here would turn
        // a geometry failure into an unphysical scattering event.
        const Vec3 original = particles.positions[i];
        const Vec3 ext = sub(geometry.bounds_max(), geometry.bounds_min());
        const double scale = std::max({1.0, std::abs(ext[0]), std::abs(ext[1]), std::abs(ext[2])});
        const double eps = std::max(options_.push_distance_a, 1e-10 * scale);
        std::array<Vec3, 4> candidates{original, original, original, original};
        const int nearest = mesh.nearest_facet(original);
        if (nearest >= 0 && nearest < static_cast<int>(mesh.facet_normals().size())) {
            const Vec3 normal = mesh.facet_normals()[static_cast<size_t>(nearest)];
            candidates[0] = add(original, mul(normal, eps));
            candidates[1] = add(original, mul(normal, -eps));
        }
        const double speed = norm(particles.velocities[i]);
        if (speed > 0.0 && std::isfinite(speed)) {
            const Vec3 direction = mul(particles.velocities[i], 1.0 / speed);
            candidates[2] = add(original, mul(direction, eps));
            candidates[3] = add(original, mul(direction, -eps));
        }
        for (const Vec3& candidate : candidates) {
            if (!geometry.contains_point(candidate)) {
                continue;
            }
            std::tie(cp, t_hit, fct) =
                geometry.trace_boundary_intersection(candidate, particles.velocities[i]);
            if (fct >= 0 && std::isfinite(t_hit)) {
                particles.positions[i] = candidate;
                collision_cache_corrections_total_.fetch_add(1, std::memory_order_relaxed);
                break;
            }
        }
    }
    if (fct < 0 || !std::isfinite(t_hit)) {
        particles.collision_positions[i] = particles.positions[i];
        particles.collision_facets[i] = -1;
        particles.collision_conditions[i] = 'R';
        particles.collision_times[i] = std::numeric_limits<double>::infinity();
        if (particles.collision_failed[static_cast<size_t>(i)] == 0) {
            particles.collision_failed[static_cast<size_t>(i)] = static_cast<std::uint8_t>(1);
            collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
        }
        return;
    }
    particles.collision_failed[static_cast<size_t>(i)] = static_cast<std::uint8_t>(0);
    particles.collision_positions[i] = cp;
    particles.collision_facets[i] = fct;
    particles.collision_conditions[i] = geometry.facet_boundary_condition(fct);
    particles.collision_times[i] = t_hit / options_.time_step_ps;
}

void TransportOperator::throw_if_collision_cache_failed(const ParticleStorage& particles,
                                                        const char* context) {
    for (size_t i = 0; i < particles.collision_failed.size(); ++i) {
        if (particles.collision_failed[i] == 0) {
            continue;
        }
        failed_ = true;
        const Vec3& p = particles.positions[i];
        if (particles.collision_failed[i] == 2) {
            throw std::runtime_error(
                std::string("Full-matrix boundary mode mapping failed while ") + context +
                " for particle " + std::to_string(i) +
                ". An inward discrete material mode and finite signed energy weight are required; "
                "no velocity-only reflection was applied.");
        }
        throw std::runtime_error(
            std::string("Boundary intersection failed while ") + context + " for particle " +
            std::to_string(i) + " at (" + std::to_string(p[0]) + ", " + std::to_string(p[1]) +
            ", " + std::to_string(p[2]) + "). Check that the STL is closed and manifold.");
    }
}

void TransportOperator::throw_if_excessive_collisions(const ParticleStorage& particles) {
    const int i = excessive_collision_particle_.load(std::memory_order_relaxed);
    if (i < 0) return;
    failed_ = true;
    const Vec3& p = particles.positions[static_cast<size_t>(i)];
    throw std::runtime_error(
        "Particle " + std::to_string(i) + " exceeded 64 transport events in one timestep at (" +
        std::to_string(p[0]) + ", " + std::to_string(p[1]) + ", " + std::to_string(p[2]) +
        "), cached facet=" + std::to_string(particles.collision_facets[i]) +
        ". Reduce simulation.time_step or inspect narrow STL features; no recovery was applied.");
}

// 函数说明：处理边界事件（透射、吸收、周期、粗糙/镜面反射）并更新粒子态。
void TransportOperator::process_boundary_collision(
    ParticleStorage& particles, int i, BoundaryLedger& ledger, int worker, std::mt19937_64& rng) {
    const auto& geometry = geometry_;
    const int facet = particles.collision_facets[i];
    const bool require_discrete_mode = options_.boundary.require_discrete_mode;
    particles.positions[i] = particles.collision_positions[i];
    phonomc::BoundaryHit hit;
    hit.valid = facet >= 0 && facet < static_cast<int>(geometry.mesh().facet_normals().size());
    hit.condition = particles.collision_conditions[i];
    if (hit.valid) {
        hit.outward_normal = geometry.mesh().facet_normals()[facet];
        hit.has_periodic_pair = geometry.has_periodic_pair(facet);
        if (hit.has_periodic_pair) {
            hit.periodic_shift = geometry.periodic_shift_for_facet(facet);
        }
        hit.rough_table = rough_sampler_.table_index(particles.material_ids[i], facet);
        // The geometry lookup scans rough facets. Keep it confined to the
        // missing-table legacy fallback, as before the event extraction.
        if (!require_discrete_mode && hit.condition == 'R' && hit.rough_table < 0) {
            hit.roughness = geometry.roughness_for_facet(facet, 0.0);
        }
    }
    const phonomc::BoundaryState incoming{particles.modes[i],
                                          particles.material_ids[i],
                                          particles.velocities[i],
                                          particles.omega[i],
                                          particles.occupation[i]};
    const auto result = phonomc::interact_boundary(
        incoming,
        hit,
        {require_discrete_mode, options_.boundary.reference_temperature},
        *materials_.at(static_cast<size_t>(particles.material_ids[i])),
        energy_weights_.at(static_cast<size_t>(particles.material_ids[i])),
        rough_sampler_,
        rng);
    record_rough_diagnostics(result.rough);
    if (result.error != phonomc::BoundaryError::None) {
        particles.collision_failed[i] = 2;
        collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    if (result.action == phonomc::BoundaryAction::Unchanged) {
        return;
    }
    if (result.action == phonomc::BoundaryAction::Absorb) {
        ledger.record_absorption(
            worker, reservoir_sampler_.reservoir_index(facet), result.absorbed_energy_ev);
        particles.alive[static_cast<size_t>(i)] = static_cast<std::uint8_t>(0);
        return;
    }
    particles.modes[i] = result.state.mode;
    particles.velocities[i] = result.state.velocity;
    particles.omega[i] = result.state.omega;
    particles.occupation[i] = result.state.occupation;
    if (result.clear_cached_energy) {
        particles.energies[i] = 0.0;
    }
    if (result.action == phonomc::BoundaryAction::Periodic) {
        particles.positions[i] = add(particles.positions[i], result.position_shift);
    }

    // Numerical placement and bookkeeping remain transport responsibilities.
    particles.positions[i] =
        add(particles.positions[i],
            mul(particles.velocities[i],
                options_.push_distance_a / std::max(norm(particles.velocities[i]), 1e-12)));
    particles.grid_ids[i] = nearest_grid_index(geometry, particles.positions[i]);
    const int destination_material = layout_.material_at(particles.positions[i]);
    if (require_discrete_mode && destination_material != particles.material_ids[i]) {
        particles.collision_failed[i] = 2;
        collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    particles.material_ids[i] = destination_material;
}

// Reject escaped particles without relocating them or changing their physical state.
void TransportOperator::validate_containment(const ParticleStorage& particles) {
    ensure_usable();
    const auto& geometry = geometry_;
    const auto& bmin = geometry.bounds_min();
    const auto& bmax = geometry.bounds_max();
    const Vec3 ext = sub(bmax, bmin);
    const double scale = std::max({1.0, std::abs(ext[0]), std::abs(ext[1]), std::abs(ext[2])});
    const double tol = 1e-10 * scale;
    std::atomic<int> escaped{-1};
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particles.size(); ++i) {
        const Vec3& p = particles.positions[static_cast<size_t>(i)];
        bool outside = false;
        for (int axis = 0; axis < 3; ++axis) {
            outside = outside || !std::isfinite(p[axis]) || p[axis] < bmin[axis] - tol ||
                      p[axis] > bmax[axis] + tol;
        }
        if (!outside && !geometry.is_box_geometry()) {
            outside = !geometry.contains_point(p);
        }
        if (outside) {
            int expected = -1;
            escaped.compare_exchange_strong(expected, i, std::memory_order_relaxed);
        }
    }
    const int i = escaped.load(std::memory_order_relaxed);
    if (i < 0) {
        return;
    }
    const Vec3& p = particles.positions[static_cast<size_t>(i)];
    failed_ = true;
    throw std::runtime_error(
        "Particle " + std::to_string(i) + " escaped the domain at (" + std::to_string(p[0]) + ", " +
        std::to_string(p[1]) + ", " + std::to_string(p[2]) +
        "), cached facet=" + std::to_string(particles.collision_facets[i]) +
        ". Refusing to relocate or rethermalize it; inspect the geometry and time step.");
}

// 函数说明：在单个时间步内推进粒子运动并处理可能的多次边界碰撞。
void TransportOperator::advance_particle(ParticleStorage& particles,
                                         int i,
                                         double dt_remaining,
                                         BoundaryLedger& ledger,
                                         int worker,
                                         std::mt19937_64& rng) {
    const auto& geometry = geometry_;
    if (particles.alive[static_cast<size_t>(i)] == 0) {
        return;
    }
    double remaining = dt_remaining;
    for (int event = 0; remaining > 1e-14 && event < 64; ++event) {
        const double t_wall = particles.collision_times[i] * options_.time_step_ps;
        const auto crossing = layout_.next_crossing(
            particles.positions[i], particles.velocities[i], particles.material_ids[i]);
        const double t_event = std::min(t_wall, crossing.time_ps);
        if (!std::isfinite(t_event) || t_event > remaining) {
            particles.positions[i] =
                add(particles.positions[i], mul(particles.velocities[i], remaining));
            if (std::isfinite(particles.collision_times[i])) {
                particles.collision_times[i] -= remaining / options_.time_step_ps;
            }
            return;
        }

        // The outer wall wins exact ties. Clamp only the elapsed flight: tiny
        // negative intersection times can arise at a numerically placed surface.
        const bool interface_first = crossing.time_ps < t_wall;
        remaining -= std::max(0.0, t_event);
        if (interface_first) {
            particles.positions[i] =
                add(particles.positions[i], mul(particles.velocities[i], crossing.time_ps));
            process_material_interface(particles, i, crossing.destination, rng);
            if (particles.collision_failed[i] != 0) return;
            const double speed = std::max(norm(particles.velocities[i]), 1e-12);
            particles.positions[i] =
                add(particles.positions[i],
                    mul(particles.velocities[i], options_.push_distance_a / speed));
            particles.grid_ids[i] = nearest_grid_index(geometry, particles.positions[i]);
        } else {
            process_boundary_collision(particles, i, ledger, worker, rng);
            if (particles.alive[i] == 0 || particles.collision_failed[i] != 0) return;
        }

        update_collision_cache_single(particles, i);
        if (particles.collision_failed[i] != 0) return;
        // Preserve the wall-only near-zero-flight correction; applying this
        // after an interface would change interface/wall event ordering.
        if (!interface_first && particles.collision_times[i] * options_.time_step_ps < 1e-12) {
            particles.positions[i] =
                add(particles.positions[i], mul(particles.velocities[i], 1e-12));
            update_collision_cache_single(particles, i);
            if (particles.collision_failed[i] != 0) return;
        }
    }
    if (remaining > 1e-14) {
        // Report after the OpenMP region. Never discard an unfinished flight.
        int expected = -1;
        excessive_collision_particle_.compare_exchange_strong(
            expected, i, std::memory_order_relaxed);
    }
}

void TransportOperator::process_material_interface(ParticleStorage& particles,
                                                   int i,
                                                   int next_material,
                                                   std::mt19937_64& rng) {
    const phonomc::BoundaryState incoming{particles.modes[i],
                                          particles.material_ids[i],
                                          particles.velocities[i],
                                          particles.omega[i],
                                          particles.occupation[i]};
    const auto result = interface_interaction_.interact(incoming, next_material, rng);
    if (result.interface_events) {
        interface_events_total_.fetch_add(result.interface_events, std::memory_order_relaxed);
    }
    if (result.interface_transmitted) {
        interface_transmitted_total_.fetch_add(result.interface_transmitted,
                                               std::memory_order_relaxed);
    }
    if (result.interface_reflected) {
        interface_reflected_total_.fetch_add(result.interface_reflected, std::memory_order_relaxed);
    }
    if (result.error != phonomc::BoundaryError::None) {
        particles.collision_failed[i] = 2;
        collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    particles.material_ids[i] = result.state.material;
    particles.modes[i] = result.state.mode;
    particles.velocities[i] = result.state.velocity;
    particles.omega[i] = result.state.omega;
    particles.occupation[i] = result.state.occupation;
}

void TransportOperator::record_rough_diagnostics(const phonomc::RoughSampleEvents& events) {
    if (events.events) {
        rough_events_total_.fetch_add(events.events, std::memory_order_relaxed);
    }
    if (events.specular) {
        rough_specular_selected_.fetch_add(events.specular, std::memory_order_relaxed);
    }
    if (events.diffuse) {
        rough_diffuse_selected_.fetch_add(events.diffuse, std::memory_order_relaxed);
    }
    if (events.residual_window) {
        rough_residual_window_selected_.fetch_add(events.residual_window,
                                                  std::memory_order_relaxed);
    }
    if (events.residual_fallback) {
        rough_residual_window_fallback_.fetch_add(events.residual_fallback,
                                                  std::memory_order_relaxed);
    }
    if (events.missing_data) {
        rough_fallback_missing_rough_data_.fetch_add(events.missing_data,
                                                     std::memory_order_relaxed);
    }
    if (events.missing_spec_match) {
        rough_fallback_missing_spec_match_.fetch_add(events.missing_spec_match,
                                                     std::memory_order_relaxed);
    }
    if (events.outgoing_pool) {
        rough_fallback_outgoing_pool_.fetch_add(events.outgoing_pool, std::memory_order_relaxed);
    }
    if (events.global_random) {
        rough_fallback_global_random_.fetch_add(events.global_random, std::memory_order_relaxed);
    }
}

TransportDiagnostics TransportOperator::diagnostics() const {
    TransportDiagnostics d;
    d.collision_cache_corrections =
        collision_cache_corrections_total_.load(std::memory_order_relaxed);
    d.collision_cache_failures = collision_cache_failures_total_.load(std::memory_order_relaxed);
    d.interface_events = interface_events_total_.load(std::memory_order_relaxed);
    d.interface_transmitted = interface_transmitted_total_.load(std::memory_order_relaxed);
    d.interface_reflected = interface_reflected_total_.load(std::memory_order_relaxed);
    d.rough.events = rough_events_total_.load(std::memory_order_relaxed);
    d.rough.specular = rough_specular_selected_.load(std::memory_order_relaxed);
    d.rough.diffuse = rough_diffuse_selected_.load(std::memory_order_relaxed);
    d.rough.residual_window = rough_residual_window_selected_.load(std::memory_order_relaxed);
    d.rough.residual_fallback = rough_residual_window_fallback_.load(std::memory_order_relaxed);
    d.rough.missing_data = rough_fallback_missing_rough_data_.load(std::memory_order_relaxed);
    d.rough.missing_spec_match = rough_fallback_missing_spec_match_.load(std::memory_order_relaxed);
    d.rough.outgoing_pool = rough_fallback_outgoing_pool_.load(std::memory_order_relaxed);
    d.rough.global_random = rough_fallback_global_random_.load(std::memory_order_relaxed);
    return d;
}
} // namespace phonomc
