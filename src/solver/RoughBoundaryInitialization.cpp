#include "MonteCarloSolver.h"
#include "PhononMaterial.h"
#include "SimulationDomain.h"

#include <algorithm>
#include <cmath>
#include <iostream>

void MonteCarloSolver::initialize_rough_boundary_scattering(const SimulationDomain& geometry) {
    std::vector<phonomc::RoughFacetSpec> facets;
    facets.reserve(static_cast<size_t>(geometry.mesh().facet_count()));
    const bool require_discrete_mode = uses_conservative_transport(args_);
    for (int f = 0; f < geometry.mesh().facet_count(); ++f) {
        const char condition = geometry.facet_boundary_condition(f);
        const bool reflecting =
            condition != 'T' && !(condition == 'P' && geometry.has_periodic_pair(f));
        facets.push_back({require_discrete_mode ? reflecting : geometry.is_rough_facet(f),
                          geometry.mesh().facet_normals()[f],
                          std::max(0.0, geometry.roughness_for_facet(f, 0.0))});
    }
    rough_sampler_.configure(facets,
                             materials_,
                             &std::cout,
                             require_discrete_mode,
                             uses_conservative_transport(args_) ? background_temperature_ : 0.0);
}
