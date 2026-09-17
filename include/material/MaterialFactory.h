#pragma once
#include "PhononMaterial.h"

struct SimulationConfig;

namespace phonomc {
// Application adapter: resolves config paths and applies the explicit test-only
// PHONOMC_ALLOW_SYNTHETIC_MATERIAL=1 fallback. Preserves startup diagnostics.
PhononMaterial load_phonon_material(const SimulationConfig& config, int material_index);
} // namespace phonomc
