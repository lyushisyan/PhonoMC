#pragma once

#include <vector>
#include "solver/CellParticleIndex.h"

class PhononMaterial;

namespace phonomc {

class ParticleStorage;

// Read-only cell statistics from the immediately preceding temperature update.
// Every cell contains one material; counts include all stored particles.
struct CollisionGridView {
    const std::vector<int>& particle_counts;
    const std::vector<int>& material_ids;
    const std::vector<double>& temperatures;
    const CellParticleIndex* cells = nullptr;
};

struct CollisionStep {
    double time_step_ps;
    double particle_volume_a3;
    bool local_lifetime_temperature;
    double lifetime_temperature_kelvin;
};

// Finite-step energy-conserving RTA. Owns reusable scratch, but no geometry,
// RNG, output stream or solver configuration. Only occupations are changed;
// the caller refreshes energies/temperatures afterwards. Not concurrently
// reentrant: one instance per solver. Internal cell loops may use OpenMP.
class RtaCollisionOperator {
public:
    double apply(
        ParticleStorage& particles,
        const std::vector<const PhononMaterial*>& materials,
        const CollisionGridView& grid,
        const CollisionStep& step);

private:
    CellParticleIndex fallback_cells_;
    std::vector<double> weights_;
};

}  // namespace phonomc
