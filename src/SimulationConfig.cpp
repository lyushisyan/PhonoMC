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

const char* to_string(InitialTemperatureMode mode) {
    return mode == InitialTemperatureMode::Linear ? "linear" : "uniform";
}

const char* to_string(BoundaryCondition condition) {
    switch (condition) {
    case BoundaryCondition::ThermalReservoir: return "T";
    case BoundaryCondition::Periodic: return "P";
    case BoundaryCondition::Rough: return "R";
    }
    return "R";
}

const char* to_string(HeatSourceProfile profile) {
    return profile == HeatSourceProfile::Gaussian ? "gaussian" : "uniform";
}

const char* to_string(HeatSourceTimeProfile profile) {
    return profile == HeatSourceTimeProfile::Square ? "square" : "constant";
}

char boundary_condition_code(BoundaryCondition condition) {
    return to_string(condition)[0];
}

namespace {
// 函数说明：清理输入文本两端空白，保证配置与数据解析的稳健性。
std::string trim(const std::string& s) {
    const auto first = std::find_if_not(s.begin(), s.end(), [](unsigned char ch) { return std::isspace(ch); });
    if (first == s.end()) {
        return "";
    }
    const auto last = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char ch) { return std::isspace(ch); }).base();
    return std::string(first, last);
}

// 函数说明：将关键字统一为小写，避免输入大小写差异影响流程。
std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return s;
}

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
    const double value = std::stod(t, &consumed);
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
    } catch (const std::exception&) {
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
std::vector<std::array<double, 6>> parse_region_array(const std::string& value, const std::string& field_name) {
    std::vector<std::array<double, 6>> boxes;
    for (const auto& row : parse_array_elements(value)) {
        const auto coords = parse_number_array(row);
        if (coords.size() != 6) {
            throw std::runtime_error(
                field_name +
                " row must be [xmin,ymin,zmin,xmax,ymax,zmax]. Point selectors [x,y,z] are no longer supported.");
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
void parse_toml_assignments(const std::string& path, std::unordered_map<std::string, std::string>& kv) {
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
            static const std::unordered_set<std::string> allowed_sections {
                "geometry", "simulation", "boundary", "heat_source", "io"
            };
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
            throw std::runtime_error("Top-level TOML keys are not supported; place '" + key + "' in its documented section.");
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
        if (!kv.emplace(full_key, value).second) {
            throw std::runtime_error("Duplicate TOML key: " + full_key);
        }
    }
}

std::optional<std::string> get_value(
    const std::unordered_map<std::string, std::string>& kv,
    const char* key) {
    const auto it = kv.find(key);
    if (it != kv.end()) {
        return it->second;
    }
    return std::nullopt;
}

void reject_unknown_keys(const std::unordered_map<std::string, std::string>& kv) {
    static const std::unordered_set<std::string> allowed {
        "geometry.model", "geometry.sizes", "geometry.merge_coplanar_facets",
        "simulation.particle_count", "simulation.time_step", "simulation.iterations",
        "simulation.convergence_write_interval", "simulation.random_seed", "simulation.compute_kappa",
        "simulation.profile_timers", "simulation.progress_temperature_summary_only",
        "simulation.temperature_lookup_dt", "simulation.background_temperature",
        "simulation.lifetime_temperature", "simulation.initial_temperature",
        "simulation.grid_xyz",
        "boundary.boundary_conditions", "boundary.boundary_values",
        "boundary.boundary_position", "boundary.periodic_pair",
        "heat_source.enabled", "heat_source.min", "heat_source.max",
        "heat_source.power_density", "heat_source.power_densities",
        "heat_source.profile", "heat_source.center", "heat_source.centers",
        "heat_source.sigma", "heat_source.half_width", "heat_source.time_profile", "heat_source.start_time",
        "heat_source.end_time", "heat_source.period", "heat_source.on_duration",
        "heat_source.duty_cycle", "heat_source.amplitude",
        "io.material_folder", "io.output_folder"
    };
    for (const auto& [key, _] : kv) {
        if (allowed.find(key) == allowed.end()) {
            throw std::runtime_error("Unknown TOML key: " + key);
        }
    }
}

// 函数说明：将区域集合扁平化为 token，供边界赋值阶段按区域匹配 facet。
std::vector<std::string> flatten_boxes(const std::string& mode, const std::vector<std::array<double, 6>>& boxes) {
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
    if (std::abs(args.heat_source_power_density) <= 0.0) {
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
    try {
        const double temperature = parse_double_scalar(raw);
        if (temperature >= 0.0) {
            return temperature;
        }
    } catch (...) {
    }
    throw std::runtime_error(
        "background_temperature must be an unquoted finite non-negative number; "
        "the local background mode has been removed.");
}

// 函数说明：寿命查询温度可使用局部温度或一个固定数值。
std::pair<bool, double> parse_lifetime_temperature_reference(
    const std::string& value,
    const std::string& field_name) {
    const std::string raw = trim(value);
    const std::string token = to_lower(unquote(raw));
    if (is_quoted(raw)) {
        if (token == "local") {
            return {true, 300.0};
        }
        throw std::runtime_error(field_name + " string value must be 'local'. Numeric fixed temperatures must not be quoted.");
    }
    try {
        size_t consumed = 0;
        const std::string number = normalise_number_token(raw);
        const double temperature = std::stod(number, &consumed);
        if (consumed == number.size() && std::isfinite(temperature) && temperature >= 0.0) {
            return {false, temperature};
        }
    } catch (...) {
    }
    throw std::runtime_error(field_name + " must be a finite non-negative number or the string 'local'.");
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
    throw std::runtime_error(
        "Unsupported boundary condition at index " + std::to_string(index) +
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
    try {
        size_t consumed = 0;
        const std::string number = normalise_number_token(token);
        const double temperature = std::stod(number, &consumed);
        if (consumed == number.size() && std::isfinite(temperature)) {
            return {InitialTemperatureMode::Uniform, temperature};
        }
    } catch (...) {
    }
    throw std::runtime_error("initial_temperature must be a finite numeric value or 'linear'.");
}

// 函数说明：校验温度参考配置，避免非法固定温度进入热物性查询。
void validate_temperature_reference_config(const SimulationConfig& args) {
    if (!std::isfinite(args.background_temperature) || args.background_temperature < 0.0) {
        throw std::runtime_error("background_temperature must be a finite non-negative value.");
    }
    if (!std::isfinite(args.lifetime_temperature) || args.lifetime_temperature < 0.0) {
        throw std::runtime_error("lifetime_temperature must be a finite non-negative value.");
    }
}

size_t flattened_region_count(const std::vector<std::string>& values, const std::string& field_name) {
    if (values.empty()) {
        return 0;
    }
    if (values.front() != "relative_box" || (values.size() - 1) % 6 != 0) {
        throw std::runtime_error(field_name + " has an invalid internal region representation.");
    }
    return (values.size() - 1) / 6;
}

void validate_boundary_config(const SimulationConfig& args) {
    const size_t selector_count = flattened_region_count(args.boundary_position, "boundary.boundary_position");
    const size_t periodic_region_count = flattened_region_count(args.periodic_pair, "boundary.periodic_pair");
    if (selector_count == 0) {
        if (!args.boundary_conditions.empty() || !args.boundary_values.empty() || periodic_region_count != 0) {
            throw std::runtime_error(
                "boundary.boundary_position is required when boundary conditions, values, or periodic pairs are provided.");
        }
        return;
    }
    if (args.boundary_conditions.size() != selector_count) {
        throw std::runtime_error(
            "boundary.boundary_conditions must contain exactly one entry per boundary_position region.");
    }
    if (args.boundary_values.size() != selector_count) {
        throw std::runtime_error(
            "boundary.boundary_values must contain exactly one entry per boundary_position region; use 0 for periodic boundaries.");
    }

    size_t periodic_condition_count = 0;
    for (size_t i = 0; i < selector_count; ++i) {
        const BoundaryCondition condition = args.boundary_conditions[i];
        const double value = args.boundary_values[i];
        if (condition == BoundaryCondition::ThermalReservoir && value < 0.0) {
            throw std::runtime_error("Thermal reservoir boundary values must be non-negative temperatures.");
        }
        if (condition == BoundaryCondition::Rough && value < 0.0) {
            throw std::runtime_error("Rough boundary values must be non-negative roughness values.");
        }
        if (condition == BoundaryCondition::Periodic) {
            ++periodic_condition_count;
            if (value != 0.0) {
                throw std::runtime_error("Periodic boundary values must be 0.");
            }
        }
    }
    if (periodic_region_count != periodic_condition_count || periodic_region_count % 2 != 0) {
        throw std::runtime_error(
            "boundary.periodic_pair must contain every periodic boundary region exactly once, in pairs.");
    }
}

void validate_heat_source_config(const SimulationConfig& args) {
    if (!args.heat_source_enabled) {
        return;
    }
    if (!(args.heat_source_power_density > 0.0) || !std::isfinite(args.heat_source_power_density)) {
        throw std::runtime_error("enabled heat_source requires a finite positive power_density.");
    }
    for (double q : args.heat_source_power_densities) {
        if (!(q > 0.0) || !std::isfinite(q)) {
            throw std::runtime_error("heat_source.power_densities values must be finite and positive.");
        }
    }
    if (args.heat_source_profile == HeatSourceProfile::Uniform) {
        if (args.heat_source_min.size() != 3 || args.heat_source_max.size() != 3) {
            throw std::runtime_error("uniform heat_source requires min and max arrays with 3 values.");
        }
    } else {
        const bool has_single_center = args.heat_source_center.size() == 3;
        const bool has_multi_centers = !args.heat_source_centers.empty();
        if ((!has_single_center && !has_multi_centers) || args.heat_source_sigma.size() != 3) {
            throw std::runtime_error("gaussian heat_source requires center=[x,y,z] or centers=[[x,y,z],...] and sigma=[sx,sy,sz].");
        }
        if (has_multi_centers && (args.heat_source_centers.size() % 3 != 0)) {
            throw std::runtime_error("heat_source.centers must contain one or more [x,y,z] rows.");
        }
        const size_t nsrc = has_multi_centers ? args.heat_source_centers.size() / 3 : 1;
        if (!args.heat_source_power_densities.empty() &&
            args.heat_source_power_densities.size() != nsrc) {
            throw std::runtime_error("heat_source.power_densities must match the number of heat_source.centers.");
        }
        if (!args.heat_source_half_width.empty() && args.heat_source_half_width.size() != 3) {
            throw std::runtime_error("heat_source.half_width must contain 3 values if provided.");
        }
        if ((!args.heat_source_min.empty() || !args.heat_source_max.empty()) &&
            (args.heat_source_min.size() != 3 || args.heat_source_max.size() != 3)) {
            throw std::runtime_error("gaussian heat_source min/max clipping box must be omitted or contain 3 values each.");
        }
    }
    if (!std::isfinite(args.heat_source_amplitude) || args.heat_source_amplitude < 0.0) {
        throw std::runtime_error("heat_source.amplitude must be finite and non-negative.");
    }
    if (!std::isfinite(args.heat_source_time_start)) {
        throw std::runtime_error("heat_source.start_time must be finite.");
    }
    if (!std::isfinite(args.heat_source_time_end)) {
        throw std::runtime_error("heat_source.end_time must be finite.");
    }
    if (!std::isfinite(args.heat_source_duty_cycle) ||
        args.heat_source_duty_cycle < 0.0 || args.heat_source_duty_cycle > 1.0) {
        throw std::runtime_error("heat_source.duty_cycle must be in [0, 1].");
    }
    if (args.heat_source_time_profile == HeatSourceTimeProfile::Square) {
        if (!(args.heat_source_period > 0.0) || !std::isfinite(args.heat_source_period)) {
            throw std::runtime_error("square heat_source.time_profile requires a finite positive period.");
        }
        if (args.heat_source_on_duration >= 0.0 &&
            (!std::isfinite(args.heat_source_on_duration) ||
             args.heat_source_on_duration > args.heat_source_period)) {
            throw std::runtime_error("heat_source.on_duration must be finite and no larger than period.");
        }
    }
}

// 函数说明：解析 TOML 输入并构造仿真配置对象（含单位与默认规则）。
SimulationConfig parse_toml_file(const std::string& path) {
    SimulationConfig args;
    args.input_directory = std::filesystem::absolute(std::filesystem::path(path)).parent_path().string();

    std::unordered_map<std::string, std::string> kv;
    parse_toml_assignments(path, kv);
    reject_unknown_keys(kv);

    if (auto v = get_value(kv, "geometry.model"); v.has_value()) {
        args.model = unquote(*v);
        if (to_lower(args.model) == "cylinder") {
            throw std::runtime_error("geometry.model='cylinder' is no longer supported; use 'box' or a mesh file.");
        }
    }
    if (auto v = get_value(kv, "geometry.sizes"); v.has_value()) {
        args.sizes = parse_number_array(*v);
        // Input geometry sizes are in nm; convert to internal Angstrom.
        for (double& x : args.sizes) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(kv, "simulation.particle_count"); v.has_value()) {
        args.particle_count = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.time_step"); v.has_value()) {
        args.time_step = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.iterations"); v.has_value()) {
        args.iterations = parse_int_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.convergence_write_interval"); v.has_value()) {
        args.convergence_write_interval = parse_int_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.random_seed"); v.has_value()) {
        args.random_seed = parse_uint64_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.compute_kappa"); v.has_value()) {
        args.compute_kappa = parse_bool_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.profile_timers"); v.has_value()) {
        args.profile_timers = parse_bool_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.progress_temperature_summary_only"); v.has_value()) {
        args.progress_temperature_summary_only = parse_bool_scalar(*v);
    }
    if (auto v = get_value(kv, "geometry.merge_coplanar_facets"); v.has_value()) {
        args.merge_coplanar_facets = parse_bool_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.temperature_lookup_dt"); v.has_value()) {
        args.temperature_lookup_dt = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "simulation.background_temperature"); v.has_value()) {
        args.background_temperature = parse_background_temperature_reference(*v);
    }
    if (auto v = get_value(kv, "simulation.lifetime_temperature"); v.has_value()) {
        const auto [is_local, temperature] =
            parse_lifetime_temperature_reference(*v, "lifetime_temperature");
        args.lifetime_temperature_is_local = is_local;
        args.lifetime_temperature = temperature;
    }

    if (auto v = get_value(kv, "simulation.initial_temperature"); v.has_value()) {
        const std::string t = trim(*v);
        if (!t.empty() && t.front() == '[') {
            throw std::runtime_error("simulation.initial_temperature must be a numeric scalar or 'linear', not an array.");
        }
        args.initial_temperature = parse_initial_temperature(t);
    }

    if (const auto gv = get_value(kv, "simulation.grid_xyz"); gv.has_value()) {
        const auto gd = parse_number_array(*gv);
        if (gd.size() != 3) {
            throw std::runtime_error("grid_xyz must contain 3 values [nx, ny, nz].");
        }
        for (double value : gd) {
            if (std::abs(value - std::round(value)) > 1e-12) {
                throw std::runtime_error("grid_xyz values must be integers.");
            }
        }
        args.grid = {
            static_cast<int>(std::llround(gd[0])),
            static_cast<int>(std::llround(gd[1])),
            static_cast<int>(std::llround(gd[2]))
        };
    } else {
        throw std::runtime_error("Missing required key: simulation.grid_xyz.");
    }

    if (auto v = get_value(kv, "boundary.boundary_conditions"); v.has_value()) {
        args.boundary_conditions = parse_boundary_conditions(parse_string_array(*v));
    }
    if (auto v = get_value(kv, "boundary.boundary_values"); v.has_value()) {
        args.boundary_values = parse_number_array(*v);
    }

    if (auto pts = get_value(kv, "boundary.boundary_position"); pts.has_value()) {
        args.boundary_position = flatten_boxes("relative_box", parse_region_array(*pts, "boundary_position"));
    }

    if (auto pts = get_value(kv, "boundary.periodic_pair"); pts.has_value()) {
        args.periodic_pair = flatten_boxes("relative_box", parse_region_array(*pts, "periodic_pair"));
    }

    if (auto v = get_value(kv, "io.material_folder"); v.has_value()) {
        args.material_folder = unquote(*v);
    }
    if (auto v = get_value(kv, "io.output_folder"); v.has_value()) {
        args.output_folder = unquote(*v);
    }

    bool heat_source_enabled_set = false;
    if (auto v = get_value(kv, "heat_source.enabled"); v.has_value()) {
        args.heat_source_enabled = parse_bool_scalar(*v);
        heat_source_enabled_set = true;
    }
    if (auto v = get_value(kv, "heat_source.min"); v.has_value()) {
        args.heat_source_min = parse_number_array(*v);
        // Heat-source min/max coordinates are specified in nm, matching
        // geometry.sizes.  Internally, PhonoMC stores geometry coordinates in
        // Angstrom.
        for (double& x : args.heat_source_min) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(kv, "heat_source.max"); v.has_value()) {
        args.heat_source_max = parse_number_array(*v);
        for (double& x : args.heat_source_max) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(kv, "heat_source.power_density"); v.has_value()) {
        args.heat_source_power_density = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "heat_source.power_densities"); v.has_value()) {
        args.heat_source_power_densities = parse_number_array(*v);
    }
    if (auto v = get_value(kv, "heat_source.profile"); v.has_value()) {
        args.heat_source_profile = parse_heat_source_profile(*v);
    }
    if (auto v = get_value(kv, "heat_source.center"); v.has_value()) {
        args.heat_source_center = parse_number_array(*v);
        // Heat-source center coordinates are specified in nm, matching geometry.sizes.
        // Internally, PhonoMC stores geometry coordinates in Angstrom.
        for (double& x : args.heat_source_center) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(kv, "heat_source.centers"); v.has_value()) {
        args.heat_source_centers = parse_vec3_array_flat(*v, "heat_source.centers");
        for (double& x : args.heat_source_centers) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(kv, "heat_source.sigma"); v.has_value()) {
        args.heat_source_sigma = parse_number_array(*v);
        // Heat-source Gaussian widths are specified in nm.  Non-positive values
        // remain non-positive after conversion and keep the "uniform axis" meaning.
        for (double& x : args.heat_source_sigma) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(kv, "heat_source.half_width"); v.has_value()) {
        args.heat_source_half_width = parse_number_array(*v);
        for (double& x : args.heat_source_half_width) {
            x *= 10.0;
        }
    }
    if (auto v = get_value(kv, "heat_source.time_profile"); v.has_value()) {
        args.heat_source_time_profile = parse_heat_source_time_profile(*v);
    }
    if (auto v = get_value(kv, "heat_source.start_time"); v.has_value()) {
        args.heat_source_time_start = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "heat_source.end_time"); v.has_value()) {
        args.heat_source_time_end = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "heat_source.period"); v.has_value()) {
        args.heat_source_period = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "heat_source.on_duration"); v.has_value()) {
        args.heat_source_on_duration = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "heat_source.duty_cycle"); v.has_value()) {
        args.heat_source_duty_cycle = parse_double_scalar(*v);
    }
    if (auto v = get_value(kv, "heat_source.amplitude"); v.has_value()) {
        args.heat_source_amplitude = parse_double_scalar(*v);
    }
    if (!heat_source_enabled_set && infer_heat_source_enabled(args)) {
        args.heat_source_enabled = true;
    }

    return args;
}

}  // namespace

// 函数说明：加载规范的分节 TOML 输入并完成统一校验。
SimulationConfig load_simulation_config(const std::string& path) {
    if (to_lower(std::filesystem::path(path).extension().string()) != ".toml") {
        throw std::runtime_error("Only sectioned TOML input files (*.toml) are supported.");
    }
    SimulationConfig args = parse_toml_file(path);
    if (!(args.temperature_lookup_dt > 0.0) || !std::isfinite(args.temperature_lookup_dt)) {
        throw std::runtime_error("temperature_lookup_dt must be a finite positive value.");
    }
    validate_temperature_reference_config(args);
    validate_boundary_config(args);
    validate_heat_source_config(args);
    if (args.convergence_write_interval <= 0) {
        throw std::runtime_error("convergence_write_interval must be a positive integer.");
    }
    if (!(args.particle_count > 0.0) || !std::isfinite(args.particle_count)) {
        throw std::runtime_error("particle_count must be a finite positive value.");
    }
    if (!(args.time_step > 0.0) || !std::isfinite(args.time_step)) {
        throw std::runtime_error("time_step must be a finite positive value.");
    }
    if (args.iterations <= 0) {
        throw std::runtime_error("iterations must be a positive integer.");
    }
    if (args.grid.nx <= 0 || args.grid.ny <= 0 || args.grid.nz <= 0) {
        throw std::runtime_error("grid_xyz values must be positive integers.");
    }
    if (args.model == "box") {
        if (args.sizes.size() != 3 ||
            std::any_of(args.sizes.begin(), args.sizes.end(), [](double value) { return !(value > 0.0); })) {
            throw std::runtime_error("geometry.sizes for box must contain three positive values.");
        }
    }
    return args;
}

// 函数说明：按序号创建结果目录，避免覆盖历史仿真输出。
std::string create_indexed_output_folder(const std::string& output_folder) {
    if (output_folder.empty()) {
        throw std::runtime_error("output_folder cannot be empty.");
    }

    namespace fs = std::filesystem;
    fs::path base_abs = fs::absolute(fs::path(output_folder)).lexically_normal();
    if (base_abs.filename().empty()) {
        base_abs = base_abs.parent_path();
    }
    const fs::path parent = base_abs.parent_path();
    const auto stem = base_abs.filename().string();

    int max_idx = -1;
    if (fs::exists(parent) && fs::is_directory(parent)) {
        for (const auto& ent : fs::directory_iterator(parent)) {
            if (!ent.is_directory()) {
                continue;
            }
            const auto name = ent.path().filename().string();
            const auto prefix = stem + "_";
            if (name.rfind(prefix, 0) != 0) {
                continue;
            }
            const auto suffix = name.substr(prefix.size());
            if (suffix.empty()) {
                continue;
            }
            bool all_digits = std::all_of(suffix.begin(), suffix.end(), [](unsigned char ch) { return std::isdigit(ch); });
            if (!all_digits) {
                continue;
            }
            max_idx = std::max(max_idx, std::stoi(suffix));
        }
    }

    const int next_idx = max_idx + 1;
    const fs::path out = parent / (stem + "_" + std::to_string(next_idx));
    fs::create_directories(out);
    return out.string();
}
