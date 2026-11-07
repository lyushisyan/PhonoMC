#pragma once
#include <random>
#include <utility>
#include <vector>

class SimulationDomain;
class PhononMaterial;

namespace phonomc {
class ReservoirSampler;
class ParticleStorage;
class BoundaryLedger;

enum class ReservoirGeneration { OneToOne, FixedIncomingFlux };

struct ReservoirInjectionOptions {
    ReservoirGeneration generation = ReservoirGeneration::OneToOne;
    // Independent policies: multimaterial flux can still use nonlinear Bose
    // occupations and the legacy velocity reversal; conservative transport cannot.
    bool conservative = false;
    bool linearized = false;
    double reference_temperature = 300;
    double push_distance_a = 1e-10;
};

// Borrows immutable geometry, materials and a configured sampler. No solver
// back-reference, owned RNG, absorption history or transport/collision state.
// Borrowed services must outlive the injector and must not be reconfigured.
class ReservoirInjector {
public:
    ReservoirInjector(ReservoirInjectionOptions options,
                      const SimulationDomain& geometry,
                      const std::vector<const PhononMaterial*>& materials,
                      std::vector<double> energy_weights,
                      const ReservoirSampler& sampler);

    // Call outside parallel regions while the ledger is collecting. OneToOne
    // reads the previous step's counts (initial rounded flux on the first step).
    // For each reservoir: count draw, all mode/time draws, all surface draws,
    // then append and record energy in that order. Returns storage index/time_ps.
    // The caller builds collision caches and advances the remaining flight time.
    // An exception may leave appended particles and ledger entries: abort the step.
    std::vector<std::pair<int, double>>
    inject(ParticleStorage& particles, BoundaryLedger& ledger, std::mt19937_64& rng) const;

private:
    ReservoirInjectionOptions options_;
    const SimulationDomain& geometry_;
    std::vector<const PhononMaterial*> materials_;
    std::vector<double> energy_weights_;
    const ReservoirSampler& sampler_;
};
} // namespace phonomc
