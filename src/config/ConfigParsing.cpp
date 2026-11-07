#include "ConfigInternal.h"
#include "ConfigText.h"
#include "SimulationConfig.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace phonomc::config_detail {
namespace {
// 函数说明：移除配置行内注释且保留引号内容，避免误删有效字段。
std::string strip_inline_comment(const std::string& line) {
    bool in_quote = false;
    char quote_char = '\0';
    bool escape = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (escape) {
            escape = false;
            continue;
        }
        if (in_quote && c == '\\') {
            escape = true;
            continue;
        }
        if ((c == '"' || c == '\'') && (!in_quote || c == quote_char)) {
            if (!in_quote) {
                in_quote = true;
                quote_char = c;
            } else {
                in_quote = false;
                quote_char = '\0';
            }
            continue;
        }
        if (!in_quote && c == '#') {
            return line.substr(0, i);
        }
    }
    return line;
}

// 函数说明：判断字符串是否被引号包裹，用于安全去引号。
bool is_quoted(const std::string& s) {
    if (s.size() < 2) {
        return false;
    }
    return (s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\'');
}

// 函数说明：去除外层引号并保留原始值语义。
std::string unquote(const std::string& s) {
    const std::string t = trim(s);
    if (is_quoted(t)) {
        return t.substr(1, t.size() - 2);
    }
    return t;
}

// 函数说明：标准化数字 token（去空白、去下划线）便于数值转换。
std::string normalise_number_token(std::string s) {
    s = trim(unquote(s));
    s.erase(std::remove(s.begin(), s.end(), '_'), s.end());
    return s;
}

// 函数说明：解析标量浮点参数，用于粒子数、时间步等基础输入。
double parse_double_scalar(const std::string& s) {
    const std::string t = normalise_number_token(s);
    if (t.empty()) {
        throw std::runtime_error("Expected numeric value, got empty token.");
    }
    size_t consumed = 0;
    double value;
    try {
        value = std::stod(t, &consumed);
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("Expected a numeric value, got: " + s);
    } catch (const std::out_of_range&) {
        throw std::runtime_error("Numeric value out of range: " + s);
    }
    if (consumed != t.size() || !std::isfinite(value)) {
        throw std::runtime_error("Expected a finite numeric value, got: " + s);
    }
    return value;
}

// 函数说明：解析严格整数参数，拒绝小数和越界值。
int parse_int_scalar(const std::string& s) {
    const double value = parse_double_scalar(s);
    const double rounded = std::round(value);
    if (std::abs(value - rounded) > 1e-12 ||
        rounded < static_cast<double>(std::numeric_limits<int>::min()) ||
        rounded > static_cast<double>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Expected an integer value, got: " + s);
    }
    return static_cast<int>(rounded);
}

// 函数说明：解析可重现随机流使用的 64 位无符号种子。
std::uint64_t parse_uint64_scalar(const std::string& s) {
    const std::string t = normalise_number_token(s);
    if (t.empty() || t.front() == '-') {
        throw std::runtime_error("random_seed must be an unsigned 64-bit integer.");
    }
    size_t consumed = 0;
    try {
        const unsigned long long value = std::stoull(t, &consumed, 10);
        if (consumed != t.size()) {
            throw std::runtime_error("random_seed must be an unsigned 64-bit integer.");
        }
        return static_cast<std::uint64_t>(value);
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("random_seed must be an unsigned 64-bit integer.");
    } catch (const std::out_of_range&) {
        throw std::runtime_error("random_seed must be an unsigned 64-bit integer.");
    }
}

// 函数说明：解析标准 TOML 布尔值。
bool parse_bool_scalar(const std::string& s) {
    const std::string t = to_lower(trim(s));
    if (t == "true") {
        return true;
    }
    if (t == "false") {
        return false;
    }
    throw std::runtime_error("Expected TOML boolean true or false, got: " + s);
}

int parse_axis(const std::string& s, const std::string& field_name) {
    const std::string value = to_lower(unquote(s));
    if (value == "x") {
        return 0;
    }
    if (value == "y") {
        return 1;
    }
    if (value == "z") {
        return 2;
    }
    throw std::runtime_error(field_name + " must be one of 'x', 'y', or 'z'.");
}

InterfaceModel parse_interface_model(const std::string& s) {
    const std::string value = to_lower(unquote(s));
    if (value == "diffuse_mismatch" || value == "dmm") {
        return InterfaceModel::DiffuseMismatch;
    }
    if (value == "acoustic_mismatch" || value == "amm") {
        return InterfaceModel::AcousticMismatch;
    }
    if (value == "mixed_mismatch" || value == "mmm") {
        return InterfaceModel::MixedMismatch;
    }
    throw std::runtime_error("materials.interface_model must be dmm, amm or mmm.");
}

// 函数说明：按顶层分隔符切分数组文本，正确处理嵌套括号与引号。
std::vector<std::string> split_top_level(const std::string& s, char delim) {
    std::vector<std::string> out;
    std::string cur;
    int depth = 0;
    bool in_quote = false;
    char quote_char = '\0';
    bool escape = false;
    for (char c : s) {
        if (escape) {
            cur.push_back(c);
            escape = false;
            continue;
        }
        if (in_quote && c == '\\') {
            cur.push_back(c);
            escape = true;
            continue;
        }
        if ((c == '"' || c == '\'') && (!in_quote || c == quote_char)) {
            if (!in_quote) {
                in_quote = true;
                quote_char = c;
            } else {
                in_quote = false;
                quote_char = '\0';
            }
            cur.push_back(c);
            continue;
        }
        if (!in_quote) {
            if (c == '[') {
                ++depth;
            } else if (c == ']') {
                --depth;
            } else if (c == delim && depth == 0) {
                out.push_back(trim(cur));
                cur.clear();
                continue;
            }
        }
        cur.push_back(c);
    }
    if (!trim(cur).empty()) {
        out.push_back(trim(cur));
    }
    return out;
}

// 函数说明：抽取数组元素文本，为后续数值/字符串/点阵解析提供输入。
std::vector<std::string> parse_array_elements(const std::string& value) {
    const std::string t = trim(value);
    if (t.size() < 2 || t.front() != '[' || t.back() != ']') {
        throw std::runtime_error("Expected TOML array value, got: " + t);
    }
    const std::string inner = trim(t.substr(1, t.size() - 2));
    if (inner.empty()) {
        return {};
    }
    return split_top_level(inner, ',');
}

// 函数说明：解析数值数组参数（如尺寸、边界值、热源范围）。
std::vector<double> parse_number_array(const std::string& value) {
    std::vector<double> out;
    for (const auto& elem : parse_array_elements(value)) {
        out.push_back(parse_double_scalar(elem));
    }
    return out;
}

// 函数说明：解析字符串数组参数（如边界类型、初始化模式）。
std::vector<std::string> parse_string_array(const std::string& value) {
    std::vector<std::string> out;
    for (const auto& elem : parse_array_elements(value)) {
        out.push_back(unquote(elem));
    }
    return out;
}

// 函数说明：解析区域数组参数（边界选择器/周期配对），每行必须为六元组。
std::vector<std::array<double, 6>> parse_region_array(const std::string& value,
                                                      const std::string& field_name) {
    std::vector<std::array<double, 6>> boxes;
    for (const auto& row : parse_array_elements(value)) {
        const auto coords = parse_number_array(row);
        if (coords.size() != 6) {
            throw std::runtime_error(field_name +
                                     " row must be [xmin,ymin,zmin,xmax,ymax,zmax]. Point "
                                     "selectors [x,y,z] are no longer supported.");
        }
        boxes.push_back({coords[0], coords[1], coords[2], coords[3], coords[4], coords[5]});
    }
    return boxes;
}

// 函数说明：解析三维点阵数组参数，每行必须为三元组。
std::vector<double> parse_vec3_array_flat(const std::string& value, const std::string& field_name) {
    std::vector<double> out;
    for (const auto& row : parse_array_elements(value)) {
        const auto coords = parse_number_array(row);
        if (coords.size() != 3) {
            throw std::runtime_error(field_name + " row must be [x,y,z].");
        }
        out.insert(out.end(), coords.begin(), coords.end());
    }
    return out;
}

// 函数说明：将数值格式化为紧凑字符串，便于输出文件可读性。
std::string format_number(double x) {
    std::ostringstream oss;
    oss << std::setprecision(17) << x;
    return oss.str();
}

// 函数说明：统计括号平衡，用于多行 TOML 值拼接的结束判定。
int bracket_balance(const std::string& s) {
    int bal = 0;
    bool in_quote = false;
    char quote_char = '\0';
    bool escape = false;
    for (char c : s) {
        if (escape) {
            escape = false;
            continue;
        }
        if (in_quote && c == '\\') {
            escape = true;
            continue;
        }
        if ((c == '"' || c == '\'') && (!in_quote || c == quote_char)) {
            if (!in_quote) {
                in_quote = true;
                quote_char = c;
            } else {
                in_quote = false;
                quote_char = '\0';
            }
            continue;
        }
        if (!in_quote) {
            if (c == '[') {
                ++bal;
            } else if (c == ']') {
                --bal;
            }
        }
    }
    return bal;
}

// 函数说明：读取分节 TOML 键值对，不接受顶层键或静默忽略的旧格式。
void parse_toml_assignments(const std::string& path,
                            std::unordered_map<std::string, std::string>& assignments) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open TOML input file: " + path);
    }
    std::string section;
    std::string line;
    while (std::getline(in, line)) {
        std::string clean = trim(strip_inline_comment(line));
        if (clean.empty()) {
            continue;
        }
        if (clean.front() == '[' && clean.back() == ']') {
            section = to_lower(trim(clean.substr(1, clean.size() - 2)));
            static const std::unordered_set<std::string> allowed_sections{"sampling",
                                                                          "convergence",
                                                                          "geometry",
                                                                          "simulation",
                                                                          "boundary",
                                                                          "heat_source",
                                                                          "mode_excitation",
                                                                          "materials",
                                                                          "io",
                                                                          "scattering",
                                                                          "driving"};
            if (allowed_sections.find(section) == allowed_sections.end()) {
                throw std::runtime_error("Unknown TOML section: [" + section + "]");
            }
            continue;
        }

        const auto eq = clean.find('=');
        if (eq == std::string::npos) {
            throw std::runtime_error("Invalid TOML assignment: " + clean);
        }
        const std::string key = to_lower(trim(clean.substr(0, eq)));
        std::string value = trim(clean.substr(eq + 1));
        if (value.empty()) {
            throw std::runtime_error("Empty TOML value for key: " + key);
        }
        if (section.empty()) {
            throw std::runtime_error("Top-level TOML keys are not supported; place '" + key +
                                     "' in its documented section.");
        }

        if (value.front() == '[') {
            while (bracket_balance(value) > 0) {
                std::string next_line;
                if (!std::getline(in, next_line)) {
                    break;
                }
                next_line = trim(strip_inline_comment(next_line));
                if (next_line.empty()) {
                    continue;
                }
                value += " " + next_line;
            }
            if (bracket_balance(value) != 0) {
                throw std::runtime_error("Unclosed TOML array for key: " + section + "." + key);
            }
        }

        const std::string full_key = section + "." + key;
        if (!assignments.emplace(full_key, value).second) {
            throw std::runtime_error("Duplicate TOML key: " + full_key);
        }
    }
}

std::optional<std::string>
get_value(const std::unordered_map<std::string, std::string>& assignments, const char* key) {
    const auto it = assignments.find(key);
    if (it != assignments.end()) {
        return it->second;
    }
    return std::nullopt;
}

void reject_unknown_keys(const std::unordered_map<std::string, std::string>& assignments) {
    static const std::unordered_set<std::string> allowed{
        "geometry.model",
        "geometry.sizes",
        "geometry.merge_coplanar_facets",
        "geometry.extruded_z",
        "simulation.particle_count",
        "simulation.time_step",
        "simulation.iterations",
        "simulation.convergence_write_interval",
        "simulation.random_seed",
        "simulation.compute_kappa",
        "simulation.write_cell_heat_flux",
        "simulation.profile_timers",
        "simulation.progress_temperature_summary_only",
        "simulation.temperature_lookup_dt",
        "simulation.background_temperature",
        "simulation.lifetime_temperature",
        "simulation.initial_temperature",
        "simulation.grid_xyz",
        "simulation.transport_axis",
        "driving.temperature_gradient",
        "sampling.resample_interval",
        "sampling.per_mode_sign",
        "driving.warm_start_ps",
        "convergence.enabled",
        "convergence.min_time_ps",
        "convergence.window_ps",
        "convergence.windows",
        "convergence.relative_tolerance",
        "convergence.absolute_tolerance",
        "scattering.model",
        "scattering.matrix_files",
        "scattering.callaway_formulation",
        "simulation.conservative_boundaries",
        "boundary.boundary_conditions",
        "boundary.boundary_values",
        "boundary.boundary_position",
        "boundary.periodic_pair",
        "heat_source.enabled",
        "heat_source.min",
        "heat_source.max",
        "heat_source.power_density",
        "heat_source.power_densities",
        "heat_source.total_power",
        "heat_source.spectrum",
        "heat_source.branches",
        "heat_source.branch_weights",
        "heat_source.frequency_min",
        "heat_source.frequency_max",
        "heat_source.frequency_center",
        "heat_source.frequency_sigma",
        "heat_source.profile",
        "heat_source.center",
        "heat_source.centers",
        "heat_source.sigma",
        "heat_source.half_width",
        "heat_source.time_profile",
        "heat_source.start_time",
        "heat_source.end_time",
        "heat_source.period",
        "heat_source.on_duration",
        "heat_source.duty_cycle",
        "heat_source.amplitude",
        "mode_excitation.enabled",
        "mode_excitation.frequency_min",
        "mode_excitation.frequency_max",
        "mode_excitation.occupation_multiplier",
        "materials.folders",
        "materials.interface_axis",
        "materials.interface_positions",
        "materials.interface_model",
        "materials.interface_frequency_bin_thz",
        "materials.interface_amm_fraction",
        "materials.interface_parallel_bin_inv_a",
        "materials.mass_densities_kg_m3",
        "io.material_folder",
        "io.output_folder"};
    for (const auto& [key, _] : assignments) {
        if (allowed.find(key) == allowed.end()) {
            throw std::runtime_error("Unknown TOML key: " + key);
        }
    }
}

// 函数说明：将区域集合扁平化为 token，供边界赋值阶段按区域匹配 facet。
std::vector<std::string> flatten_boxes(const std::string& mode,
                                       const std::vector<std::array<double, 6>>& boxes) {
    std::vector<std::string> out;
    out.reserve(1 + boxes.size() * 6);
    out.push_back(mode);
    for (const auto& b : boxes) {
        for (int i = 0; i < 6; ++i) {
            out.push_back(format_number(b[static_cast<size_t>(i)]));
        }
    }
    return out;
}

// 函数说明：在未显式给出 enabled 时，根据 profile 与关键参数推断是否启用热源。
bool infer_heat_source_enabled(const SimulationConfig& args) {
    if (!args.heat_source_total_power.has_value() &&
        std::abs(args.heat_source_power_density) <= 0.0) {
        return false;
    }
    if (args.heat_source_profile == HeatSourceProfile::Gaussian) {
        return (args.heat_source_center.size() == 3 || args.heat_source_centers.size() >= 3) &&
               args.heat_source_sigma.size() == 3;
    }
    return args.heat_source_min.size() == 3 && args.heat_source_max.size() == 3;
}

// 函数说明：背景温度是全局不变的能量参考，只接受非负数值。
double parse_background_temperature_reference(const std::string& value) {
    const std::string raw = trim(value);
    if (is_quoted(raw)) {
        throw std::runtime_error(
            "background_temperature must be an unquoted finite non-negative number; "
            "the local background mode has been removed.");
    }
    return parse_double_scalar(raw);
}

// 函数说明：寿命查询温度可使用局部温度或一个固定数值。
std::pair<bool, double> parse_lifetime_temperature_reference(const std::string& value,
                                                             const std::string& field_name) {
    const std::string raw = trim(value);
    const std::string token = to_lower(unquote(raw));
    if (is_quoted(raw)) {
        if (token == "local") {
            return {true, 300.0};
        }
        throw std::runtime_error(
            field_name +
            " string value must be 'local'. Numeric fixed temperatures must not be quoted.");
    }
    return {false, parse_double_scalar(raw)};
}

BoundaryCondition parse_boundary_condition(const std::string& value, size_t index) {
    const std::string token = to_lower(trim(unquote(value)));
    if (token == "t") {
        return BoundaryCondition::ThermalReservoir;
    }
    if (token == "p") {
        return BoundaryCondition::Periodic;
    }
    if (token == "r") {
        return BoundaryCondition::Rough;
    }
    throw std::runtime_error("Unsupported boundary condition at index " + std::to_string(index) +
                             ": '" + value + "'. Only 'T', 'P', 'R' are supported.");
}

std::vector<BoundaryCondition> parse_boundary_conditions(const std::vector<std::string>& values) {
    std::vector<BoundaryCondition> result;
    result.reserve(values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        result.push_back(parse_boundary_condition(values[i], i));
    }
    return result;
}

HeatSourceProfile parse_heat_source_profile(const std::string& value) {
    const std::string profile = to_lower(trim(unquote(value)));
    if (profile == "uniform") {
        return HeatSourceProfile::Uniform;
    }
    if (profile == "gaussian") {
        return HeatSourceProfile::Gaussian;
    }
    throw std::runtime_error("heat_source.profile must be either 'uniform' or 'gaussian'.");
}

HeatSourceTimeProfile parse_heat_source_time_profile(const std::string& value) {
    const std::string profile = to_lower(trim(unquote(value)));
    if (profile == "constant") {
        return HeatSourceTimeProfile::Constant;
    }
    if (profile == "square") {
        return HeatSourceTimeProfile::Square;
    }
    throw std::runtime_error("heat_source.time_profile must be either 'constant' or 'square'.");
}

InitialTemperatureConfig parse_initial_temperature(const std::string& value) {
    const std::string raw = trim(value);
    const std::string token = to_lower(unquote(raw));
    if (token == "linear") {
        return {InitialTemperatureMode::Linear, 300.0};
    }
    if (is_quoted(raw)) {
        throw std::runtime_error("Numeric initial_temperature values must not be quoted.");
    }
    return {InitialTemperatureMode::Uniform, parse_double_scalar(token)};
}

void parse_material_settings(const std::unordered_map<std::string, std::string>& assignments,
                             SimulationConfig& args) {
    const auto material_folders_value = get_value(assignments, "materials.folders");
    const bool has_other_material_keys =
        get_value(assignments, "materials.interface_axis").has_value() ||
        get_value(assignments, "materials.interface_positions").has_value() ||
        get_value(assignments, "materials.interface_model").has_value() ||
        get_value(assignments, "materials.interface_frequency_bin_thz").has_value() ||
        get_value(assignments, "materials.interface_amm_fraction").has_value() ||
        get_value(assignments, "materials.interface_parallel_bin_inv_a").has_value() ||
        get_value(assignments, "materials.mass_densities_kg_m3").has_value();
    if (!material_folders_value.has_value() && has_other_material_keys) {
        throw std::runtime_error(
            "materials.folders is required when any other [materials] key is present.");
    }
    if (material_folders_value.has_value()) {
        args.material_folders = parse_string_array(*material_folders_value);
        if (args.material_folders.empty()) {
            throw std::runtime_error("materials.folders must not be empty.");
        }
    }
    if (auto v = get_value(assignments, "materials.interface_axis"); v.has_value()) {
        args.material_interface_axis = parse_axis(*v, "materials.interface_axis");
    }
    if (auto v = get_value(assignments, "materials.interface_positions"); v.has_value()) {
        args.material_interface_positions = parse_number_array(*v);
        for (double& position : args.material_interface_positions) {
            position *= 10.0;
        }
    }
    if (auto v = get_value(assignments, "materials.interface_model"); v.has_value()) {
        args.material_interface_model = parse_interface_model(*v);
    }
    if (auto v = get_value(assignments, "materials.interface_frequency_bin_thz"); v.has_value()) {
        args.material_interface_frequency_bin_thz = parse_double_scalar(*v);
    }

    if (auto v = get_value(assignments, "materials.interface_amm_fraction"); v.has_value()) {
        args.material_interface_amm_fraction = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "materials.interface_parallel_bin_inv_a"); v.has_value()) {
        args.material_interface_parallel_bin_inv_a = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "materials.mass_densities_kg_m3"); v.has_value()) {
        args.material_mass_densities_kg_m3 = parse_number_array(*v);
    }
    if (args.material_interface_model == InterfaceModel::MixedMismatch &&
        !get_value(assignments, "materials.interface_amm_fraction").has_value()) {
        throw std::runtime_error("MMM requires explicit materials.interface_amm_fraction.");
    }
    if (args.material_interface_model == InterfaceModel::AcousticMismatch &&
        get_value(assignments, "materials.interface_amm_fraction").has_value() &&
        args.material_interface_amm_fraction != 1) {
        throw std::runtime_error(
            "AMM explicit fraction must equal 1; use MMM for intermediate weights.");
    }
}

void parse_heat_source_settings(const std::unordered_map<std::string, std::string>& assignments,
                                SimulationConfig& args) {
    bool heat_source_enabled_set = false;
    if (auto v = get_value(assignments, "heat_source.enabled"); v.has_value()) {
        args.heat_source_enabled = parse_bool_scalar(*v);
        heat_source_enabled_set = true;
    }
    if (auto v = get_value(assignments, "heat_source.min"); v.has_value()) {
        args.heat_source_min = parse_number_array(*v);
    }
    if (auto v = get_value(assignments, "heat_source.max"); v.has_value()) {
        args.heat_source_max = parse_number_array(*v);
    }
    if (auto v = get_value(assignments, "heat_source.total_power"); v.has_value()) {
        if (get_value(assignments, "heat_source.power_density").has_value() ||
            get_value(assignments, "heat_source.power_densities").has_value()) {
            throw std::runtime_error("Specify either total_power or density inputs, not both.");
        }
        args.heat_source_total_power = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.spectrum"); v.has_value()) {
        const auto value = to_lower(unquote(*v));
        if (value == "thermal") {
            args.heat_source_spectrum = HeatSourceSpectrum::Thermal;
        } else if (value == "weighted") {
            args.heat_source_spectrum = HeatSourceSpectrum::Weighted;
        } else if (value == "gaussian") {
            args.heat_source_spectrum = HeatSourceSpectrum::Gaussian;
        } else {
            throw std::runtime_error(
                "heat_source.spectrum must be 'thermal', 'weighted', or 'gaussian'.");
        }
    }
    if (auto v = get_value(assignments, "heat_source.branches"); v.has_value()) {
        for (double branch : parse_number_array(*v)) {
            if (!std::isfinite(branch) || branch < 0 || branch != std::floor(branch) ||
                branch > static_cast<double>(std::numeric_limits<int>::max())) {
                throw std::runtime_error(
                    "heat_source.branches must contain nonnegative integer indices.");
            }
            args.heat_source_branches.push_back(static_cast<int>(branch));
        }
    }
    if (auto v = get_value(assignments, "heat_source.branch_weights"); v.has_value()) {
        args.heat_source_branch_weights = parse_number_array(*v);
    }
    if (auto v = get_value(assignments, "heat_source.power_density"); v.has_value()) {
        args.heat_source_power_density = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.power_densities"); v.has_value()) {
        args.heat_source_power_densities = parse_number_array(*v);
    }
    if (auto v = get_value(assignments, "heat_source.frequency_center")) {
        args.heat_source_frequency_center = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.frequency_sigma")) {
        args.heat_source_frequency_sigma = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.frequency_min"); v.has_value()) {
        args.heat_source_frequency_min = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.frequency_max"); v.has_value()) {
        args.heat_source_frequency_max = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.profile"); v.has_value()) {
        args.heat_source_profile = parse_heat_source_profile(*v);
    }
    if (auto v = get_value(assignments, "heat_source.center"); v.has_value()) {
        args.heat_source_center = parse_number_array(*v);
    }
    if (auto v = get_value(assignments, "heat_source.centers"); v.has_value()) {
        args.heat_source_centers = parse_vec3_array_flat(*v, "heat_source.centers");
    }
    if (auto v = get_value(assignments, "heat_source.sigma"); v.has_value()) {
        args.heat_source_sigma = parse_number_array(*v);
    }
    if (auto v = get_value(assignments, "heat_source.half_width"); v.has_value()) {
        args.heat_source_half_width = parse_number_array(*v);
    }
    if (auto v = get_value(assignments, "heat_source.time_profile"); v.has_value()) {
        args.heat_source_time_profile = parse_heat_source_time_profile(*v);
    }
    if (auto v = get_value(assignments, "heat_source.start_time"); v.has_value()) {
        args.heat_source_time_start = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.end_time"); v.has_value()) {
        args.heat_source_time_end = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.period"); v.has_value()) {
        args.heat_source_period = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.on_duration"); v.has_value()) {
        args.heat_source_on_duration = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.duty_cycle"); v.has_value()) {
        args.heat_source_duty_cycle = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "heat_source.amplitude"); v.has_value()) {
        args.heat_source_amplitude = parse_double_scalar(*v);
    }
    if (!heat_source_enabled_set && infer_heat_source_enabled(args)) {
        args.heat_source_enabled = true;
    }
    // All source lengths enter as nm and leave as Angstrom. Non-positive
    // Gaussian widths retain their documented uniform-axis meaning.
    for (auto* lengths : {&args.heat_source_min, &args.heat_source_max,
                          &args.heat_source_center, &args.heat_source_centers,
                          &args.heat_source_sigma, &args.heat_source_half_width}) {
        for (double& length : *lengths) length *= 10.0;
    }
}

} // namespace

// 函数说明：解析 TOML 输入并构造仿真配置对象（含单位与默认规则）。
SimulationConfig parse_toml_file(const std::string& path) {
    SimulationConfig args;
    args.input_directory =
        std::filesystem::absolute(std::filesystem::path(path)).parent_path().string();

    std::unordered_map<std::string, std::string> assignments;
    parse_toml_assignments(path, assignments);
    reject_unknown_keys(assignments);

    if (auto v = get_value(assignments, "geometry.model"); v.has_value()) {
        args.model = unquote(*v);
        if (to_lower(args.model) == "cylinder") {
            throw std::runtime_error(
                "geometry.model='cylinder' is no longer supported; use 'box' or a mesh file.");
        }
    }
    if (auto v = get_value(assignments, "geometry.sizes"); v.has_value()) {
        args.sizes = parse_number_array(*v);
        // Input geometry sizes are in nm; convert to internal Angstrom.
        for (double& x : args.sizes) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(assignments, "simulation.particle_count"); v.has_value()) {
        args.particle_count = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.time_step"); v.has_value()) {
        args.time_step = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.iterations"); v.has_value()) {
        args.iterations = parse_int_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.convergence_write_interval"); v.has_value()) {
        args.convergence_write_interval = parse_int_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.random_seed"); v.has_value()) {
        args.random_seed = parse_uint64_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.compute_kappa"); v.has_value()) {
        args.compute_kappa = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.transport_axis"); v.has_value()) {
        const auto axis = to_lower(unquote(trim(*v)));
        if (axis != "x" && axis != "y" && axis != "z") {
            throw std::runtime_error("simulation.transport_axis must be x, y, or z.");
        }
        args.transport_axis = axis == "x" ? 0 : axis == "y" ? 1 : 2;
    }
    if (auto v = get_value(assignments, "driving.temperature_gradient"); v.has_value()) {
        args.temperature_gradient = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "sampling.resample_interval")) {
        args.resample_interval = parse_int_scalar(*v);
    }
    if (auto v = get_value(assignments, "sampling.per_mode_sign")) {
        args.resample_per_mode_sign = parse_int_scalar(*v);
    }
    if (auto v = get_value(assignments, "driving.warm_start_ps")) {
        args.gradient_warm_start_ps = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "convergence.enabled")) {
        args.convergence_stop = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "convergence.min_time_ps")) {
        args.convergence_min_time_ps = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "convergence.window_ps")) {
        args.convergence_window_ps = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "convergence.windows")) {
        args.convergence_windows = parse_int_scalar(*v);
    }
    if (auto v = get_value(assignments, "convergence.relative_tolerance")) {
        args.convergence_relative_tolerance = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "convergence.absolute_tolerance")) {
        args.convergence_absolute_tolerance = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.profile_timers"); v.has_value()) {
        args.profile_timers = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.progress_temperature_summary_only");
        v.has_value()) {
        args.progress_temperature_summary_only = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "geometry.merge_coplanar_facets"); v.has_value()) {
        args.merge_coplanar_facets = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.temperature_lookup_dt"); v.has_value()) {
        args.temperature_lookup_dt = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.background_temperature"); v.has_value()) {
        args.background_temperature = parse_background_temperature_reference(*v);
    }
    if (auto v = get_value(assignments, "simulation.lifetime_temperature"); v.has_value()) {
        const auto [is_local, temperature] =
            parse_lifetime_temperature_reference(*v, "lifetime_temperature");
        args.lifetime_temperature_is_local = is_local;
        args.lifetime_temperature = temperature;
    }
    if (auto v = get_value(assignments, "scattering.model"); v.has_value()) {
        const auto model = to_lower(unquote(trim(*v)));
        if (model == "rta") {
            args.collision_model = CollisionModel::Rta;
        } else if (model == "callaway") {
            args.collision_model = CollisionModel::Callaway;
        } else {
            throw std::runtime_error("scattering.model must be 'rta' or 'callaway'.");
        }
    }
    if (auto v = get_value(assignments, "scattering.callaway_formulation")) {
        const auto formulation = to_lower(unquote(trim(*v)));
        if (formulation != "linearized" && formulation != "nonlinear") {
            throw std::runtime_error("callaway_formulation must be linearized or nonlinear");
        }
        args.callaway_nonlinear = formulation == "nonlinear";
    }
    if (auto v = get_value(assignments, "geometry.extruded_z")) {
        args.extruded_z = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.write_cell_heat_flux")) {
        args.write_cell_heat_flux = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "simulation.conservative_boundaries")) {
        args.conservative_boundaries = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "scattering.matrix_files"); v.has_value()) {
        args.scattering_matrix_files = parse_string_array(*v);
    }

    if (auto v = get_value(assignments, "simulation.initial_temperature"); v.has_value()) {
        const std::string t = trim(*v);
        if (!t.empty() && t.front() == '[') {
            throw std::runtime_error("simulation.initial_temperature must be a numeric scalar or "
                                     "'linear', not an array.");
        }
        args.initial_temperature = parse_initial_temperature(t);
    }

    if (const auto gv = get_value(assignments, "simulation.grid_xyz"); gv.has_value()) {
        const auto gd = parse_number_array(*gv);
        if (gd.size() != 3) {
            throw std::runtime_error("grid_xyz must contain 3 values [nx, ny, nz].");
        }
        for (double value : gd) {
            if (std::abs(value - std::round(value)) > 1e-12 || value < 1 ||
                value > std::numeric_limits<int>::max()) {
                throw std::runtime_error("grid_xyz values must be positive representable integers.");
            }
        }
        args.grid = {static_cast<int>(std::llround(gd[0])),
                     static_cast<int>(std::llround(gd[1])),
                     static_cast<int>(std::llround(gd[2]))};
    } else {
        throw std::runtime_error("Missing required key: simulation.grid_xyz.");
    }

    if (auto v = get_value(assignments, "boundary.boundary_conditions"); v.has_value()) {
        args.boundary_conditions = parse_boundary_conditions(parse_string_array(*v));
    }
    if (auto v = get_value(assignments, "boundary.boundary_values"); v.has_value()) {
        args.boundary_values = parse_number_array(*v);
    }

    if (auto pts = get_value(assignments, "boundary.boundary_position"); pts.has_value()) {
        args.boundary_position =
            flatten_boxes("relative_box", parse_region_array(*pts, "boundary_position"));
    }

    if (auto pts = get_value(assignments, "boundary.periodic_pair"); pts.has_value()) {
        args.periodic_pair =
            flatten_boxes("relative_box", parse_region_array(*pts, "periodic_pair"));
    }

    if (auto v = get_value(assignments, "io.material_folder"); v.has_value()) {
        args.material_folder = unquote(*v);
    }
    if (auto v = get_value(assignments, "io.output_folder"); v.has_value()) {
        args.output_folder = unquote(*v);
    }
    parse_material_settings(assignments, args);
    parse_heat_source_settings(assignments, args);
    if (auto v = get_value(assignments, "mode_excitation.enabled"); v.has_value()) {
        args.mode_excitation_enabled = parse_bool_scalar(*v);
    }
    if (auto v = get_value(assignments, "mode_excitation.frequency_min"); v.has_value()) {
        args.mode_excitation_frequency_min = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "mode_excitation.frequency_max"); v.has_value()) {
        args.mode_excitation_frequency_max = parse_double_scalar(*v);
    }
    if (auto v = get_value(assignments, "mode_excitation.occupation_multiplier"); v.has_value()) {
        args.mode_excitation_occupation_multiplier = parse_double_scalar(*v);
    }

    return args;
}

} // namespace phonomc::config_detail
