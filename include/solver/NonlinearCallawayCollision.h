#pragma once
#include "material/CallawayOperator.h"
#include "solver/RtaCollisionOperator.h"
namespace phonomc {
// Positive, finite-step moment-matched Bose BGK updates on existing carriers.
// N(dt/2), R(dt), N(dt/2); modal rates are frozen at the old cell temperature.
std::vector<double> nonlinear_callaway_evolve(const CallawayData& data,
    const std::vector<double>& occupation, double temperature, double dt_ps);
class NonlinearCallawayCollision {
public:
    void configure(const std::vector<const PhononMaterial*>& materials, double reference);
    double apply(ParticleStorage& particles, const std::vector<const PhononMaterial*>& materials,
                 const CollisionGridView& grid, const CollisionStep& step);
private:
    std::vector<const PhononMaterial*> materials_;
    CellParticleIndex fallback_cells_;
};
}
