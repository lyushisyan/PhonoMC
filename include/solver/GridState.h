#pragma once

#include "solver/CellParticleIndex.h"
#include "solver/ParticleStorage.h"

namespace phonomc {
// Mutable mesh observables owned by one solver. Geometry and material tables
// remain external. Cell membership is reused by statistics, RTA and sources.
struct GridState {
    std::vector<double> temperatures;
    std::vector<int> material_ids;
    std::vector<double> energy_density;
    std::vector<Vec3> heat_flux;
    CellParticleIndex cells;
};
}  // namespace phonomc
