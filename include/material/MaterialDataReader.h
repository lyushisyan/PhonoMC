#pragma once
#include "material/MaterialData.h"
#include <filesystem>
#include <ostream>

namespace phonomc {
LatticeData read_poscar_lattice(const std::filesystem::path& path);
// Explicit folder only: config resolution and synthetic fallback are caller policy.
// Throws std::runtime_error on invalid/missing input; returns only a complete value.
// Optional progress preserves the legacy material-file announcement.
MaterialData read_material_data(const std::filesystem::path& folder, std::ostream* progress = nullptr);
} // namespace phonomc
