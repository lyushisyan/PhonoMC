#pragma once
#include "solver/BoundaryOutcome.h"

class PhononMaterial;

namespace phonomc {
struct BoundaryHit {
    bool valid = true;
    char condition = 'R';
    std::array<double, 3> outward_normal{};
    bool has_periodic_pair = false;
    std::array<double, 3> periodic_shift{};
    int rough_table = -1;
    double roughness = 0;
};

struct BoundaryPolicy {
    bool require_discrete_mode = false;
    double reference_temperature = 300;
};

// Samples physics only. Caller applies position pushes, storage/cache mutations,
// absorption bookkeeping and diagnostics. RNG is consumed in legacy order.
BoundaryOutcome interact_boundary(const BoundaryState& incoming,
                                  const BoundaryHit& hit,
                                  const BoundaryPolicy& policy,
                                  const PhononMaterial& material,
                                  double energy_weight,
                                  const RoughBoundarySampler& sampler,
                                  std::mt19937_64& rng);
} // namespace phonomc
