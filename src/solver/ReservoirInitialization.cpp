#include "MonteCarloSolver.h"
#include "PhononMaterial.h"
#include "SimulationDomain.h"

#include <iostream>

// 函数说明：预计算热库注入概率与计数器，驱动边界粒子注入机制。
void MonteCarloSolver::initialize_reservoir_injection(const SimulationDomain& geometry) {
    const auto& mesh = geometry.mesh();
    std::vector<phonomc::ReservoirFacetSpec> facets;
    for (int facet : geometry.reservoir_facets()) {
        facets.push_back({facet,
                          material_index_at(mesh.facet_centroids()[facet]),
                          geometry.reservoir_value_for_facet(facet, 300.0),
                          mesh.facet_areas()[facet],
                          mesh.facet_normals()[facet]});
    }
    std::vector<double> densities;
    for (double volume : particles_.material_spatial_volumes_a3) {
        densities.push_back(1.0 / volume);
    }
    reservoir_sampler_.configure(
        facets, materials_, mesh.facet_count(), particle_density_, args_.time_step, densities);
    std::vector<int> initial_counts(static_cast<size_t>(reservoir_sampler_.size()));
    for (int r = 0; r < reservoir_sampler_.size(); ++r) {
        initial_counts[r] = reservoir_sampler_.initial_emission_count(r);
    }
    boundary_ledger_.configure(initial_counts);

    phonomc::ReservoirInjectionOptions options;
    options.generation = (uses_conservative_transport(args_) || materials_.size() > 1)
                             ? phonomc::ReservoirGeneration::FixedIncomingFlux
                             : phonomc::ReservoirGeneration::OneToOne;
    options.conservative = uses_conservative_transport(args_);
    options.linearized = uses_linearized_transport(args_);
    options.reference_temperature = background_temperature_;
    options.push_distance_a = push_eps_;
    reservoir_injector_ = std::make_unique<phonomc::ReservoirInjector>(
        options, geometry, materials_, material_energy_weights(), reservoir_sampler_);
    std::cout << "[init] reservoir_gen="
              << (options.generation == phonomc::ReservoirGeneration::FixedIncomingFlux
                      ? "fixed_incoming_flux"
                      : "one_to_one")
              << '\n';
}
