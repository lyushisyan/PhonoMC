#pragma once
#include "SimulationConfig.h" // InterfaceModel; no SimulationConfig instance is retained.
#include "solver/AcousticInterfaceSampler.h"
#include "solver/BoundaryOutcome.h"
#include "solver/InterfaceModeSampler.h"

namespace phonomc {
struct InterfaceInteractionConfig {
    int axis = 0;
    bool enabled = false;
    bool linearized = false;
    InterfaceModel model = InterfaceModel::DiffuseMismatch;
    double frequency_bin_thz = 0.05;
    double amm_fraction = 0;
    double parallel_bin_inv_a = 0.05;
    double reference_temperature = 300;
    std::vector<double> mass_densities_kg_m3;
};

// Owns immutable sampling tables and scalar policy; borrows immutable materials.
// No geometry, particle mutation, ledger, global counters or owned RNG.
class InterfaceInteraction {
public:
    void configure(const InterfaceInteractionConfig& config,
                   const std::vector<const PhononMaterial*>& materials,
                   const std::vector<double>& energy_weights);
    BoundaryOutcome
    interact(const BoundaryState& incoming, int next_material, std::mt19937_64& rng) const;

    const AcousticInterfaceSampler::Audit& audit(int face) const { return acoustic_.audit(face); }

private:
    InterfaceInteractionConfig config_;
    std::vector<const PhononMaterial*> materials_;
    std::vector<double> energy_weights_;
    InterfaceModeSampler diffuse_;
    AcousticInterfaceSampler acoustic_;
};
} // namespace phonomc
