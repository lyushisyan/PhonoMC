#pragma once
#include "solver/BoundaryInteraction.h"
#include "solver/ParticleStorage.h"

#include <atomic>
#include <random>
#include <utility>
#include <vector>

class SimulationDomain;
class PhononMaterial;

namespace phonomc {
class LayeredMaterialLayout;
class InterfaceInteraction;
class ReservoirSampler;
class BoundaryLedger;

struct TransportOptions {
    double time_step_ps = 1;
    double push_distance_a = 1e-10;
    BoundaryPolicy boundary;
};

// Borrow the original RNGs: serial builds use serial; OpenMP builds use the
// executing worker's stream, including stream zero in a one-thread team.
struct TransportRandomStreams {
    std::mt19937_64& serial;
    std::vector<std::mt19937_64>& workers;
    std::mt19937_64& for_worker(int worker) const;
    void validate() const;
};

struct TransportDiagnostics {
    long long collision_cache_corrections = 0;
    long long collision_cache_failures = 0;
    RoughSampleEvents rough;
    long long interface_events = 0;
    long long interface_transmitted = 0;
    long long interface_reflected = 0;
};

// Borrows immutable geometry/material/layout/sampling services for its lifetime.
// Owns transport diagnostics and failure state, but no particles, RNG or ledger.
// Public operations run serially outside parallel regions and never append or compact
// particles. The caller invalidates cell membership before moving particles.
// A failed operation may have advanced part of the population: do not resume it.
class TransportOperator {
public:
    TransportOperator(TransportOptions options,
                      const SimulationDomain& geometry,
                      const LayeredMaterialLayout& layout,
                      const std::vector<const PhononMaterial*>& materials,
                      std::vector<double> energy_weights,
                      const RoughBoundarySampler& rough,
                      const InterfaceInteraction& interfaces,
                      const ReservoirSampler& reservoirs);

    void begin_step();
    void update_collision_cache(ParticleStorage& particles, const std::vector<int>& indices);
    void advance_existing(ParticleStorage& particles,
                          BoundaryLedger& ledger,
                          TransportRandomStreams rngs);
    // Each pair is a live storage index and its emission time within this step.
    // Indices must remain valid until this batch ends; refresh their caches first.
    void advance_injected(ParticleStorage& particles,
                          const std::vector<std::pair<int, double>>& injected,
                          BoundaryLedger& ledger,
                          TransportRandomStreams rngs);
    void validate_containment(const ParticleStorage& particles);
    TransportDiagnostics diagnostics() const;

    bool failed() const { return failed_; }

private:
    void ensure_usable() const;
    void finish_advance(const ParticleStorage& particles, const char* context);
    void update_collision_cache_single(ParticleStorage& particles, int i);
    void throw_if_collision_cache_failed(const ParticleStorage& particles,
                                         const char* context);
    void throw_if_excessive_collisions(const ParticleStorage& particles);
    void advance_particle(ParticleStorage& particles,
                          int i,
                          double remaining,
                          BoundaryLedger& ledger,
                          int worker,
                          std::mt19937_64& rng);
    void process_boundary_collision(ParticleStorage& particles,
                                    int i,
                                    BoundaryLedger& ledger,
                                    int worker,
                                    std::mt19937_64& rng);
    void process_material_interface(ParticleStorage& particles,
                                    int i,
                                    int next_material,
                                    std::mt19937_64& rng);
    void record_rough_diagnostics(const RoughSampleEvents& events);

    TransportOptions options_;
    const SimulationDomain& geometry_;
    const LayeredMaterialLayout& layout_;
    std::vector<const PhononMaterial*> materials_;
    std::vector<double> energy_weights_;
    const RoughBoundarySampler& rough_sampler_;
    const InterfaceInteraction& interface_interaction_;
    const ReservoirSampler& reservoir_sampler_;
    bool failed_ = false;

    std::atomic<long long> rough_events_total_{0};
    std::atomic<long long> rough_specular_selected_{0};
    std::atomic<long long> rough_diffuse_selected_{0};
    std::atomic<long long> rough_residual_window_selected_{0};
    std::atomic<long long> rough_residual_window_fallback_{0};
    std::atomic<long long> rough_fallback_missing_rough_data_{0};
    std::atomic<long long> rough_fallback_missing_spec_match_{0};
    std::atomic<long long> rough_fallback_outgoing_pool_{0};
    std::atomic<long long> rough_fallback_global_random_{0};
    std::atomic<long long> interface_events_total_{0};
    std::atomic<long long> interface_transmitted_total_{0};
    std::atomic<long long> interface_reflected_total_{0};
    std::atomic<long long> collision_cache_corrections_total_{0};
    std::atomic<long long> collision_cache_failures_total_{0};
    std::atomic<int> excessive_collision_particle_{-1};
};
} // namespace phonomc
