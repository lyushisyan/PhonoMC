#pragma once
#include "material/CallawayOperator.h"
#include "solver/RtaCollisionOperator.h"
namespace phonomc {
// Collision quadrature is the EXISTING carrier list, with one equal statistical
// weight per carrier. Duplicate modes are distinct quadrature samples.
// Updates occupations only; never creates, removes, moves or changes a mode.
class CarrierCollisionOperator {
public:
    void configure(const std::vector<const PhononMaterial*>& materials,
                   double reference_temperature, bool rta = false);
    double apply(ParticleStorage& particles, const std::vector<const PhononMaterial*>& materials,
                 const CollisionGridView& grid, double dt_ps, double particle_volume_a3);
private:
    std::vector<const PhononMaterial*> materials_;
    std::vector<CallawayData> tables_;
    std::vector<std::vector<double>> background_;
    CellParticleIndex fallback_cells_;
};
}
