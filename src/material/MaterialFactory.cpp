#include "material/MaterialFactory.h"
#include "material/MaterialDataReader.h"
#include "SimulationConfig.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace phonomc {
namespace {
std::filesystem::path resolve_folder(const SimulationConfig& args, int mat_index) {
    namespace fs = std::filesystem;
    fs::path folder;
    if (!args.material_folders.empty()) {
        if (mat_index < 0 || mat_index >= static_cast<int>(args.material_folders.size()))
            throw std::runtime_error("Material index is outside materials.folders.");
        folder = args.material_folders[static_cast<size_t>(mat_index)];
    } else if (!args.material_folder.empty()) {
        folder = args.material_folder;
    } else {
        folder = "Material/Si";
    }
    if (!folder.is_absolute()) {
        const fs::path p1 = fs::current_path() / folder;
        const fs::path p2 = fs::path(args.input_directory) / folder;
        folder = fs::exists(p1) ? p1 : p2;
    }
    return folder.lexically_normal();
}
} // namespace

PhononMaterial load_phonon_material(const SimulationConfig& args, int mat_index) {
    // The legacy config API defaults invalid spacing to 0.1 K; the new direct
    // in-memory API instead rejects invalid explicit construction parameters.
    const double dt = std::isfinite(args.temperature_lookup_dt) && args.temperature_lookup_dt > 0
        ? args.temperature_lookup_dt : 0.1;
    std::string folder_path;
    auto material = [&]() -> PhononMaterial {
        try {
            const auto folder = resolve_folder(args, mat_index);
            folder_path = folder.string();
            std::cout << "PhononMaterial material folder: " << folder_path << '\n';
            return PhononMaterial(read_material_data(folder, &std::cout), dt,
                                  uses_conservative_transport(args) || args.heat_source_spectrum==HeatSourceSpectrum::Gaussian);
        } catch (const std::exception& ex) {
            const char* allow_fallback = std::getenv("PHONOMC_ALLOW_SYNTHETIC_MATERIAL");
            if (allow_fallback == nullptr || std::string(allow_fallback) != "1")
                throw std::runtime_error(
                    std::string("Failed to load HDF5 phonon data (") + ex.what() + "). "
                    "Simulation is aborted to avoid invalid thermal conductivity results. "
                    "Fix material_folder / HDF5 / POSCAR paths, or set PHONOMC_ALLOW_SYNTHETIC_MATERIAL=1 for test-only fallback.");
            std::cerr << "Warning: failed to load HDF5 phonon data (" << ex.what()
                      << "). PHONOMC_ALLOW_SYNTHETIC_MATERIAL=1 is set, falling back to synthetic mode bank.\n";
            if (!folder_path.empty())
                std::cout << "PhononMaterial material folder: " << folder_path << '\n';
            else if (!args.material_folder.empty())
                std::cout << "PhononMaterial material folder: " << args.material_folder << '\n';
            return PhononMaterial::synthetic_for_testing(dt);
        }
    }();
    std::cout << "PhononMaterial initialized: active_mode_count=" << material.active_mode_count() << '\n';
    return material;
}
} // namespace phonomc

// Compatibility entry point; the file-free model library never calls this adapter.
PhononMaterial::PhononMaterial(const SimulationConfig& args, int mat_index)
    : PhononMaterial(phonomc::load_phonon_material(args, mat_index)) {}
