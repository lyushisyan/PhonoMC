#pragma once

#include "material/ScatteringMatrix.h"
#include "material/CallawayOperator.h"
#include "solver/RtaCollisionOperator.h"

#include <random>
#include <vector>

class PhononMaterial;

namespace phonomc {
class ParticleStorage;

struct FullMatrixCollisionResult {
    double residual_ev = 0.0;
    std::vector<int> appended_indices;
};

// Cellwise full linearized collision operator. The matrix acts on SUMS of
// signed carrier energies per mode. Existing carriers share each modal
// increment; an initially unrepresented target mode receives a new carrier.
// Callers refresh the geometry cache for appended_indices and grid statistics.
// No geometry ownership or concurrent reentrancy; signed statistical
// occupations are permitted and must also be supported by boundary/source code.
class FullMatrixCollisionOperator {
public:
    void configure(std::vector<ScatteringMatrix> matrices,
        const std::vector<const PhononMaterial*>& materials, double reference_temperature);
    // Reuse the modal gather/scatter path with a matrix-free Callaway kernel.
    void configure_callaway(const std::vector<const PhononMaterial*>& materials,
                            double reference_temperature);
    void configure_linearized_rta(const std::vector<const PhononMaterial*>& materials,
                                 double reference_temperature);

    // resolve_carrier_decay additionally attenuates within-mode carrier
    // variation using the relaxation diagonal. Used by periodic gradient runs;
    // the default retains the legacy cell-modal increment convention.
    FullMatrixCollisionResult apply(ParticleStorage& particles,
        const std::vector<const PhononMaterial*>& materials,
        const CollisionGridView& grid, double time_step_ps, double particle_volume_a3,
        double background_temperature, std::mt19937_64& rng, bool resolve_carrier_decay = false);

private:
    void configure_relaxation(const std::vector<const PhononMaterial*>& materials,
                              double reference_temperature, bool rta);
    struct ModalCache {
        std::vector<int> flat_to_matrix;
        std::vector<double> hw, background, relaxation_rate;
    };
    std::vector<ScatteringMatrix> matrices_;
    std::vector<CallawayOperator> callaway_;
    std::vector<const PhononMaterial*> materials_;
    std::vector<ModalCache> modal_cache_;
    double reference_temperature_ = 0.0;
    CellParticleIndex fallback_cells_;
};
} // namespace phonomc
