#pragma once
#include "solver/PhononSourceOperator.h"

class SimulationDomain;
namespace phonomc {
// Homogeneous gradient source -c_lambda v_lambda.axis G, in energy/ps.
// Evaluated on existing equal-weight carrier quadrature, with its sampled
// energy component projected out in each nonempty cell. Never creates carriers.
// Signed deviations are measured from the affine local equilibrium reference;
// only the periodic departure from that reference is transported.
class TemperatureGradientDrive {
public:
    void configure(const SimulationConfig& config, const SimulationDomain& geometry,
                   const PhononMaterial& material, double particle_volume_a3);
    SourceDeposit apply(ParticleStorage& particles, const CellParticleIndex& cells,
                        double time_ps, std::mt19937_64& rng) const;
    SourceDeposit warm_start(ParticleStorage& particles, const CellParticleIndex& cells,
                            double age_ps, std::mt19937_64& rng) const;
private:
    const PhononMaterial* material_ = nullptr;
    std::vector<double> rate_, hw_, background_, capacity_, cell_scale_;
    std::vector<std::array<double,3>> centers_;
    std::array<double,3> half_width_{};
    double reference_ = 0, energy_weight_ = 0;
};
}
