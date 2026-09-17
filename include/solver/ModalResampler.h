#pragma once
#include "solver/BackgroundCache.h"
#include "solver/CellParticleIndex.h"
#include "solver/ParticleStorage.h"
#include <random>

namespace phonomc {
struct ResampleResult { int removed = 0; double energy_residual_ev = 0; };
// Systematic resampling within (cell, physical mode, sign). No sign cancellation,
// mode changes or centroid relocation. Selected carriers retain valid ray caches.
ResampleResult resample_modal_carriers(ParticleStorage& particles, CellParticleIndex& cells,
    int cell_count, const PhononMaterial& material, const BackgroundCache& background,
    int per_sign, std::mt19937_64& rng);
}
