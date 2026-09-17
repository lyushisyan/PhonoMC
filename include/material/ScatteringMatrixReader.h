#pragma once
#include "material/ScatteringMatrix.h"
#include <filesystem>

namespace phonomc {

// Read the version-1 PhonoMC sparse energy-generator interchange schema.
// Raw phono3py collision files must first pass through the converter, which
// resolves their basis, units and full-grid mode ordering explicitly.
// Throws std::runtime_error for missing, malformed or nonconserving input.
ScatteringMatrixData read_scattering_matrix(const std::filesystem::path& path);

} // namespace phonomc
