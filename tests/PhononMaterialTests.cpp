#include "PhononMaterial.h"
#include "SimulationConfig.h"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

int main() {
    const auto tag = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const fs::path root = fs::temp_directory_path() / ("phonomc-material-tests-" + std::to_string(tag));
    fs::create_directories(root);

#if defined(_WIN32)
    _putenv_s("PHONOMC_ALLOW_SYNTHETIC_MATERIAL", "");
#else
    unsetenv("PHONOMC_ALLOW_SYNTHETIC_MATERIAL");
#endif

    SimulationConfig cfg;
    cfg.input_directory = root.string();
    cfg.material_folder = root.string();

    bool rejected = false;
    try {
        (void) PhononMaterial(cfg, 0);
    } catch (const std::runtime_error& ex) {
        rejected = std::string(ex.what()).find("POSCAR") != std::string::npos;
    }

    std::error_code ec;
    fs::remove_all(root, ec);
    if (!rejected) {
        std::cerr << "Missing POSCAR was not rejected by strict material loading\n";
        return 1;
    }
    std::cout << "Strict POSCAR loading test passed\n";
    return 0;
}
