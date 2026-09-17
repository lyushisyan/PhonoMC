#include "MonteCarloSolver.h"
#include "PhononMaterial.h"
#include "SimulationDomain.h"
#include <algorithm>
#include <cmath>
#include <iostream>

void MonteCarloSolver::initialize_rough_boundary_scattering(const SimulationDomain& geometry) {
    std::vector<phonomc::RoughFacetSpec> facets;
    facets.reserve(static_cast<size_t>(geometry.mesh().facet_count()));
    const bool full_matrix = uses_conservative_transport(args_);
    for (int f = 0; f < geometry.mesh().facet_count(); ++f) {
        const char condition = geometry.facet_boundary_condition(f);
        const bool reflecting = condition != 'T' && !(condition == 'P' && geometry.has_periodic_pair(f));
        facets.push_back({full_matrix ? reflecting : geometry.is_rough_facet(f), geometry.mesh().facet_normals()[f],
                          std::max(0.0, geometry.roughness_for_facet(f, 0.0))});
    }
    rough_sampler_.configure(facets, materials_, &std::cout, full_matrix,
        uses_conservative_transport(args_) ? background_temperature_ : 0.0);
}

std::array<int, 2> MonteCarloSolver::select_reflected_mode(
    const PhononMaterial& phonon, int rough_idx, const std::array<int, 2>& in_mode,
    double& out_occupation, double in_occupation) const {
    const auto reflected = rough_sampler_.sample(phonon, rough_idx, in_mode, in_occupation,
                                                background_temperature_reference(), thread_rng(),
                                                uses_conservative_transport(args_));
    if (reflected.diagnostics.events)
        rough_events_total_.fetch_add(reflected.diagnostics.events, std::memory_order_relaxed);
    if (reflected.diagnostics.specular)
        rough_specular_selected_.fetch_add(reflected.diagnostics.specular, std::memory_order_relaxed);
    if (reflected.diagnostics.diffuse)
        rough_diffuse_selected_.fetch_add(reflected.diagnostics.diffuse, std::memory_order_relaxed);
    if (reflected.diagnostics.residual_window)
        rough_residual_window_selected_.fetch_add(reflected.diagnostics.residual_window, std::memory_order_relaxed);
    if (reflected.diagnostics.residual_fallback)
        rough_residual_window_fallback_.fetch_add(reflected.diagnostics.residual_fallback, std::memory_order_relaxed);
    if (reflected.diagnostics.missing_data)
        rough_fallback_missing_rough_data_.fetch_add(reflected.diagnostics.missing_data, std::memory_order_relaxed);
    if (reflected.diagnostics.missing_spec_match)
        rough_fallback_missing_spec_match_.fetch_add(reflected.diagnostics.missing_spec_match, std::memory_order_relaxed);
    if (reflected.diagnostics.outgoing_pool)
        rough_fallback_outgoing_pool_.fetch_add(reflected.diagnostics.outgoing_pool, std::memory_order_relaxed);
    if (reflected.diagnostics.global_random)
        rough_fallback_global_random_.fetch_add(reflected.diagnostics.global_random, std::memory_order_relaxed);
    out_occupation = reflected.occupation;
    return reflected.mode;
}

// 函数说明：根据粗糙度与入射条件计算镜面反射概率。
double MonteCarloSolver::compute_roughness_specularity(const SimulationDomain& geometry, const PhononMaterial& phonon, int i, int facet) const {
    const double eta = std::max(0.0, geometry.roughness_for_facet(facet, 0.0));
    const double speed = std::max(1e-12, norm(particles_.velocities[i]));
    const Vec3 n = geometry.mesh().facet_normals()[facet];
    const double incidence_cos = std::min(1.0, std::max(0.0, std::abs(dot(particles_.velocities[i], n)) / speed));
    const double k_norm = std::max(1e-12, phonon.mode_angular_frequency(particles_.modes[i]) / speed);
    const double x = 2.0 * eta * incidence_cos * k_norm;
    const double p = std::exp(-(x * x));
    return std::clamp(p, 0.0, 1.0);
}
