#include "MonteCarloSolver.h"

#include "PhononMaterial.h"
#include "SimulationDomain.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// 函数说明：预计算热库注入概率与计数器，驱动边界粒子注入机制。
void MonteCarloSolver::initialize_reservoir_injection(const SimulationDomain& geometry) {
    const auto& mesh = geometry.mesh();
    std::vector<phonomc::ReservoirFacetSpec> facets;
    for (int facet : geometry.reservoir_facets()) {
        facets.push_back({facet, material_index_at(mesh.facet_centroids()[facet]),
            geometry.reservoir_value_for_facet(facet, 300.0), mesh.facet_areas()[facet],
            mesh.facet_normals()[facet]});
    }
    std::vector<double> densities;
    for (double volume:particles_.material_spatial_volumes_a3) densities.push_back(1.0/volume);
    reservoir_sampler_.configure(facets, materials_, mesh.facet_count(), particle_density_, time_step_, densities);
    std::vector<int> initial_counts(static_cast<size_t>(reservoir_sampler_.size()));
    for (int r = 0; r < reservoir_sampler_.size(); ++r)
        initial_counts[r] = reservoir_sampler_.initial_emission_count(r);
    boundary_ledger_.configure(initial_counts);
    std::cout << "[init] reservoir_gen="
              << ((uses_conservative_transport(args_) || materials_.size()>1) ? "fixed_incoming_flux" : "one_to_one") << '\n';
}

// 函数说明：从热库边界按统计规则注入新粒子并返回注入时刻偏移。
std::vector<std::pair<int, double>> MonteCarloSolver::inject_particles_from_reservoirs(const SimulationDomain& geometry) {
    std::vector<std::pair<int, double>> inserted;
    if (reservoir_sampler_.size() <= 0) {
        return inserted;
    }

    auto& rng = thread_rng();
    const auto& mesh = geometry.mesh();

    for (int r = 0; r < reservoir_sampler_.size(); ++r) {
        const auto& entry = reservoir_sampler_.facet(r);
        const int material_id = entry.material;
        const PhononMaterial& phonon = material(material_id);
        const auto& modes = phonon.active_mode_list();
        // Multimaterial DMM needs prescribed incoming phase-space density,
        // independent of the material-dependent outgoing carrier population.
        const int n_emit = (uses_conservative_transport(args_) || materials_.size()>1)
            ? reservoir_sampler_.sample_flux_emission_count(r, rng) : boundary_ledger_.emission_count(r);
        const auto emissions = reservoir_sampler_.sample(r, n_emit, rng);
        if (emissions.empty()) continue;
        // Keep position sampling after the complete mode/time batch.
        std::vector<Vec3> surf = mesh.sample_surface_points(
            static_cast<int>(emissions.size()), std::vector<int>{entry.facet}, rng);
        const Vec3 n = entry.outward_normal;
        const double Tres = entry.temperature;

        for (size_t k = 0; k < emissions.size(); ++k) {
            const std::array<int, 2> mode = modes[emissions[k].active_index];
            Vec3 gv = phonon.mode_group_velocity(mode);
            if (dot(gv, n) > 0.0) {
                if (uses_conservative_transport(args_))
                    throw std::runtime_error("Full-matrix reservoir sampled an outward mode.");
                gv = mul(gv, -1.0);
            }
            const double speed = std::max(1e-12, norm(gv));
            Vec3 pos = add(surf[k], mul(gv, push_eps_ / speed));

            const double occupation = uses_linearized_transport(args_)
                ? phonon.linearized_bose_occupation(Tres, background_temperature_, mode)
                : phonon.bose_occupation(Tres, mode);
            phonomc::ParticleRecord particle;
            particle.mode = mode;
            particle.material_id = material_id;
            particle.position = pos;
            particle.velocity = gv;
            particle.collision_position = pos;
            particle.temperature = Tres;
            particle.omega = phonon.mode_angular_frequency(mode);
            particle.occupation = occupation;
            particle.grid_id = nearest_grid_index(geometry, pos);
            const int new_idx = particles_.append(particle);
            inserted.push_back({new_idx, emissions[k].time_ps});

            constexpr double kHbarEvPs = 6.582119569e-4;
            const double neq = phonon.bose_occupation(background_temperature_reference(), mode);
            const double raw_deviation_ev = kHbarEvPs * phonon.mode_angular_frequency(mode) * (occupation - neq);
            boundary_ledger_.record_injection(raw_deviation_ev * material_energy_weight(material_id));
        }
    }
    return inserted;
}
