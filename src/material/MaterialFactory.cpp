#include "material/MaterialFactory.h"

#include "SimulationConfig.h"
#include "material/MaterialDataReader.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace phonomc {
namespace {
std::filesystem::path resolve_folder(const SimulationConfig& args, int material_index) {
    namespace fs = std::filesystem;
    fs::path folder;
    if (!args.material_folders.empty()) {
        if (material_index < 0 ||
            material_index >= static_cast<int>(args.material_folders.size())) {
            throw std::runtime_error("Material index is outside materials.folders.");
        }
        folder = args.material_folders[static_cast<size_t>(material_index)];
    } else if (!args.material_folder.empty()) {
        folder = args.material_folder;
    } else {
        folder = "Material/Si";
    }
    if (!folder.is_absolute()) {
        const fs::path working_directory_candidate = fs::current_path() / folder;
        const fs::path input_directory_candidate = fs::path(args.input_directory) / folder;
        folder = fs::exists(working_directory_candidate) ? working_directory_candidate
                                                         : input_directory_candidate;
    }
    return folder.lexically_normal();
}
} // namespace

PhononMaterial load_phonon_material(const SimulationConfig& args, int material_index) {
    // Reject invalid configuration before file I/O or the opt-in test fixture path.
    if (!std::isfinite(args.temperature_lookup_dt) || args.temperature_lookup_dt <= 0.0) {
        throw std::invalid_argument("temperature_lookup_dt must be finite and positive (K).");
    }
    const double lookup_step_k = args.temperature_lookup_dt;
    std::string folder_path;
    auto material = [&]() -> PhononMaterial {
        try {
            const auto folder = resolve_folder(args, material_index);
            folder_path = folder.string();
            std::cout << "PhononMaterial material folder: " << folder_path << '\n';
            return PhononMaterial(read_material_data(folder, &std::cout),
                                  lookup_step_k,
                                  uses_conservative_transport(args) ||
                                      args.heat_source_spectrum == HeatSourceSpectrum::Gaussian,
                                  args.collision_model == CollisionModel::Callaway);
        } catch (const std::exception& ex) {
            const char* allow_fallback = std::getenv("PHONOMC_ALLOW_SYNTHETIC_MATERIAL");
            if (allow_fallback == nullptr || std::string(allow_fallback) != "1") {
                throw std::runtime_error(
                    std::string("Failed to load HDF5 phonon data (") + ex.what() +
                    "). "
                    "Simulation is aborted to avoid invalid thermal conductivity results. "
                    "Fix material_folder / HDF5 / POSCAR paths, or set "
                    "PHONOMC_ALLOW_SYNTHETIC_MATERIAL=1 for test-only fallback.");
            }
            std::cerr << "Warning: failed to load HDF5 phonon data (" << ex.what()
                      << "). PHONOMC_ALLOW_SYNTHETIC_MATERIAL=1 is set, falling back to synthetic "
                         "mode bank.\n";
            if (!folder_path.empty()) {
                std::cout << "PhononMaterial material folder: " << folder_path << '\n';
            } else if (!args.material_folder.empty()) {
                std::cout << "PhononMaterial material folder: " << args.material_folder << '\n';
            }
            return PhononMaterial::synthetic_for_testing(lookup_step_k);
        }
    }();
    std::cout << "PhononMaterial initialized: active_mode_count=" << material.active_mode_count()
              << ", bz_refinement=" << material.bz_refinement() << '\n';
    return material;
}
} // namespace phonomc

// Compatibility entry point; the file-free model library never calls this adapter.
PhononMaterial::PhononMaterial(const SimulationConfig& args, int material_index)
    : PhononMaterial(phonomc::load_phonon_material(args, material_index)) {}
