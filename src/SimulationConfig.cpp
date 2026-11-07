#include "SimulationConfig.h"

#include "config/ConfigInternal.h"
#include "config/ConfigText.h"

#include <filesystem>
#include <stdexcept>

const char* to_string(InitialTemperatureMode mode) {
    return mode == InitialTemperatureMode::Linear ? "linear" : "uniform";
}

const char* to_string(BoundaryCondition condition) {
    switch (condition) {
    case BoundaryCondition::ThermalReservoir:
        return "T";
    case BoundaryCondition::Periodic:
        return "P";
    case BoundaryCondition::Rough:
        return "R";
    }
    return "R";
}

const char* to_string(HeatSourceProfile profile) {
    return profile == HeatSourceProfile::Gaussian ? "gaussian" : "uniform";
}

const char* to_string(HeatSourceTimeProfile profile) {
    return profile == HeatSourceTimeProfile::Square ? "square" : "constant";
}

const char* to_string(InterfaceModel model) {
    switch (model) {
    case InterfaceModel::DiffuseMismatch:
        return "diffuse_mismatch";
    case InterfaceModel::AcousticMismatch:
        return "acoustic_mismatch";
    case InterfaceModel::MixedMismatch:
        return "mixed_mismatch";
    }
    return "diffuse_mismatch";
}

char boundary_condition_code(BoundaryCondition condition) {
    return to_string(condition)[0];
}

const char* to_string(HeatSourceSpectrum spectrum) {
    if (spectrum == HeatSourceSpectrum::Gaussian) {
        return "gaussian";
    }
    return spectrum == HeatSourceSpectrum::Weighted ? "weighted" : "thermal";
}

// 函数说明：加载规范的分节 TOML 输入并完成统一校验。
SimulationConfig load_simulation_config(const std::string& path) {
    if (phonomc::config_detail::to_lower(std::filesystem::path(path).extension().string()) !=
        ".toml") {
        throw std::runtime_error("Only sectioned TOML input files (*.toml) are supported.");
    }
    SimulationConfig args = phonomc::config_detail::parse_toml_file(path);
    phonomc::config_detail::validate_loaded_config(args);
    return args;
}
