#pragma once
#include <array>
#include <filesystem>
#include <vector>

namespace phonomc {
// Coordinates remain in file units and face/vertex order is preserved.
// Domain setup owns path resolution, STL nm-to-Angstrom scaling and origin shifts.
struct SurfaceMeshData {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::array<int, 3>> faces;
};

SurfaceMeshData read_surface_mesh(const std::filesystem::path& path);
} // namespace phonomc
