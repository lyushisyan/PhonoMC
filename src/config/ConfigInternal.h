#pragma once
#include "SimulationConfig.h"

namespace phonomc::config_detail {
// Application facade: parse syntax and convert input units, then validate physics.
SimulationConfig parse_toml_file(const std::string& path);
void validate_loaded_config(const SimulationConfig& config);
} // namespace phonomc::config_detail
