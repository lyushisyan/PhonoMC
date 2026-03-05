#include "SimulationConfig.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

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

// 函数说明：将旧版参数文件拆分为 token 序列，兼容历史输入格式。
std::vector<std::string> tokenize_legacy(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Failed to open input file: " + path);
    }

    std::vector<std::string> tokens;
    std::string line;
    while (std::getline(in, line)) {
        const auto clean = trim(line);
        if (clean.empty() || clean[0] == '#') {
            continue;
        }
        std::istringstream iss(clean);
        std::string tk;
        while (iss >> tk) {
            tokens.push_back(tk);
        }
    }
    return tokens;
}

// 函数说明：识别旧版命令行风格选项前缀（--）。
bool is_option(const std::string& tk) {
    return tk.rfind("--", 0) == 0;
}

// 函数说明：根据扩展名判断输入文件是否按 TOML 语法解析。
bool looks_like_toml(const std::string& path) {
    const std::string ext = to_lower(std::filesystem::path(path).extension().string());
    return ext == ".toml";
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
    return std::stod(t);
}

// 函数说明：解析整数参数并进行四舍五入兼容。
int parse_int_scalar(const std::string& s) {
    return static_cast<int>(std::llround(parse_double_scalar(s)));
}

// 函数说明：解析布尔开关，支持 true/false 与数值兼容写法。
bool parse_bool_scalar(const std::string& s) {
    std::string t = to_lower(trim(unquote(s)));
    if (t == "true" || t == "yes" || t == "on") {
        return true;
    }
    if (t == "false" || t == "no" || t == "off") {
        return false;
    }
    if (t.empty()) {
        return false;
    }
    return std::abs(parse_double_scalar(t)) > 0.0;
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

// 函数说明：解析三维点数组参数（边界点、周期配对点）。
std::vector<std::array<double, 3>> parse_point_array(const std::string& value) {
    std::vector<std::array<double, 3>> points;
    for (const auto& row : parse_array_elements(value)) {
        const auto coords = parse_number_array(row);
        if (coords.size() != 3) {
            throw std::runtime_error("Point array must contain triplets [x,y,z].");
        }
        points.push_back({coords[0], coords[1], coords[2]});
    }
    return points;
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

// 函数说明：读取 TOML 键值对并展开段名为扁平键，统一后续取值逻辑。
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
            continue;
        }

        const auto eq = clean.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        const std::string key = to_lower(trim(clean.substr(0, eq)));
        std::string value = trim(clean.substr(eq + 1));
        if (value.empty()) {
            continue;
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
        }

        if (!section.empty()) {
            kv[section + "." + key] = value;
            if (kv.find(key) == kv.end()) {
                kv[key] = value;
            }
        } else {
            kv[key] = value;
        }
    }
}

// 函数说明：按候选键顺序从键值表中查找第一个有效值。
std::optional<std::string> get_first_value(
    const std::unordered_map<std::string, std::string>& kv,
    std::initializer_list<const char*> keys) {
    for (const char* k : keys) {
        const auto it = kv.find(to_lower(k));
        if (it != kv.end()) {
            return it->second;
        }
    }
    return std::nullopt;
}

// 函数说明：将点集合扁平化为旧接口需要的顺序 token。
std::vector<std::string> flatten_points(const std::string& mode, const std::vector<std::array<double, 3>>& pts) {
    std::vector<std::string> out;
    out.reserve(1 + pts.size() * 3);
    out.push_back(mode);
    for (const auto& p : pts) {
        out.push_back(format_number(p[0]));
        out.push_back(format_number(p[1]));
        out.push_back(format_number(p[2]));
    }
    return out;
}

// 函数说明：解析 TOML 输入并构造仿真配置对象（含单位与默认规则）。
SimulationConfig parse_toml_file(const std::string& path) {
    SimulationConfig args;
    args.input_directory = std::filesystem::absolute(std::filesystem::path(path)).parent_path().string();

    std::unordered_map<std::string, std::string> kv;
    parse_toml_assignments(path, kv);

    if (auto v = get_first_value(kv, {"geometry.model", "model"}); v.has_value()) {
        args.model = unquote(*v);
    }
    if (auto v = get_first_value(kv, {"geometry.sizes", "sizes"}); v.has_value()) {
        args.sizes = parse_number_array(*v);
        // Input geometry sizes are in nm; convert to internal Angstrom.
        for (double& x : args.sizes) {
            x *= 10.0;
        }
    }
    if (auto v = get_first_value(kv, {"simulation.particle_count", "particle_count"}); v.has_value()) {
        args.particle_count = parse_double_scalar(*v);
    }
    if (auto v = get_first_value(kv, {"simulation.time_step", "time_step"}); v.has_value()) {
        args.time_step = parse_double_scalar(*v);
    }
    if (auto v = get_first_value(kv, {"simulation.iterations", "iterations"}); v.has_value()) {
        args.iterations = parse_int_scalar(*v);
    }
    if (auto v = get_first_value(kv, {
            "simulation.compute_kappa",
            "compute_kappa"}); v.has_value()) {
        args.compute_kappa = parse_bool_scalar(*v);
    }

    if (auto v = get_first_value(kv, {"simulation.initial_temperature", "initial_temperature"}); v.has_value()) {
        const std::string t = trim(*v);
        if (!t.empty() && t.front() == '[') {
            args.initial_temperature = parse_string_array(t);
        } else {
            args.initial_temperature = {unquote(t)};
        }
    }

    if (const auto gv = get_first_value(kv, {"simulation.grid_xyz", "grid_xyz"}); gv.has_value()) {
        const auto gd = parse_number_array(*gv);
        if (gd.size() != 3) {
            throw std::runtime_error("grid_xyz must contain 3 values [nx, ny, nz].");
        }
        args.grid_layout = {
            "grid",
            std::to_string(static_cast<int>(std::llround(gd[0]))),
            std::to_string(static_cast<int>(std::llround(gd[1]))),
            std::to_string(static_cast<int>(std::llround(gd[2])))
        };
    } else {
        throw std::runtime_error("Missing required key: simulation.grid_xyz (or grid_xyz).");
    }

    if (auto v = get_first_value(kv, {"boundary.boundary_conditions", "boundary_conditions"}); v.has_value()) {
        args.boundary_conditions = parse_string_array(*v);
    }
    if (auto v = get_first_value(kv, {"boundary.boundary_values", "boundary_values"}); v.has_value()) {
        args.boundary_values = parse_number_array(*v);
    }

    if (auto pts = get_first_value(kv, {
            "boundary.boundary_position", "boundary_position"}); pts.has_value()) {
        args.boundary_position = flatten_points("relative", parse_point_array(*pts));
    }

    if (auto pts = get_first_value(kv, {
            "boundary.periodic_pair", "periodic_pair"}); pts.has_value()) {
        args.periodic_pair = flatten_points("relative", parse_point_array(*pts));
    }

    if (auto v = get_first_value(kv, {"io.material_folder", "material_folder"}); v.has_value()) {
        args.material_folder = unquote(*v);
    }
    if (auto v = get_first_value(kv, {
            "io.output_folder", "output_folder"}); v.has_value()) {
        args.output_folder = unquote(*v);
    }

    bool heat_source_enabled_set = false;
    if (auto v = get_first_value(kv, {"heat_source.enabled"}); v.has_value()) {
        args.heat_source_enabled = parse_bool_scalar(*v);
        heat_source_enabled_set = true;
    }
    if (auto v = get_first_value(kv, {"heat_source.min"}); v.has_value()) {
        args.heat_source_min = parse_number_array(*v);
    }
    if (auto v = get_first_value(kv, {"heat_source.max"}); v.has_value()) {
        args.heat_source_max = parse_number_array(*v);
    }
    if (auto v = get_first_value(kv, {"heat_source.power_density"}); v.has_value()) {
        args.heat_source_power_density = parse_double_scalar(*v);
    }
    if (!heat_source_enabled_set &&
        args.heat_source_min.size() == 3 &&
        args.heat_source_max.size() == 3 &&
        std::abs(args.heat_source_power_density) > 0.0) {
        args.heat_source_enabled = true;
    }

    return args;
}

// 函数说明：解析旧版参数文件并映射到统一配置对象。
SimulationConfig parse_legacy_file(const std::string& path) {
    SimulationConfig args;
    args.input_directory = std::filesystem::absolute(std::filesystem::path(path)).parent_path().string();
    const auto tokens = tokenize_legacy(path);

    std::unordered_map<std::string, std::vector<std::string>> kv;
    std::string current_key;

    for (const auto& tk : tokens) {
        if (is_option(tk)) {
            current_key = tk;
            kv[current_key] = {};
        } else if (!current_key.empty()) {
            kv[current_key].push_back(tk);
        }
    }

    if (auto it = kv.find("--model"); it != kv.end() && !it->second.empty()) {
        args.model = it->second.front();
    }
    if (auto it = kv.find("--sizes"); it != kv.end() && !it->second.empty()) {
        args.sizes.clear();
        for (const auto& v : it->second) {
            args.sizes.push_back(std::stod(v) * 10.0);
        }
    }
    if (auto it = kv.find("--particle_count"); it != kv.end() && !it->second.empty()) {
        args.particle_count = std::stod(it->second.front());
    }
    if (auto it = kv.find("--time_step"); it != kv.end() && !it->second.empty()) {
        args.time_step = std::stod(it->second.front());
    }
    if (auto it = kv.find("--iterations"); it != kv.end() && !it->second.empty()) {
        args.iterations = std::stoi(it->second.front());
    }
    if (auto it = kv.find("--compute_kappa"); it != kv.end() && !it->second.empty()) {
        args.compute_kappa = parse_bool_scalar(it->second.front());
    }
    if (auto it = kv.find("--grid_xyz"); it != kv.end() && it->second.size() >= 3) {
        args.grid_layout = {"grid", it->second[0], it->second[1], it->second[2]};
    } else {
        throw std::runtime_error("Missing required option: --grid_xyz nx ny nz");
    }
    if (auto it = kv.find("--initial_temperature"); it != kv.end()) {
        args.initial_temperature = it->second;
    }
    if (auto it = kv.find("--boundary_conditions"); it != kv.end()) {
        args.boundary_conditions = it->second;
    }
    if (auto it = kv.find("--boundary_position"); it != kv.end()) {
        args.boundary_position = it->second;
    }
    if (auto it = kv.find("--boundary_values"); it != kv.end()) {
        args.boundary_values.clear();
        for (const auto& v : it->second) {
            args.boundary_values.push_back(std::stod(v));
        }
    }
    if (auto it = kv.find("--periodic_pair"); it != kv.end()) {
        args.periodic_pair = it->second;
    }
    if (auto it = kv.find("--material_folder"); it != kv.end() && !it->second.empty()) {
        args.material_folder = it->second.front();
    }
    if (auto it = kv.find("--output_folder"); it != kv.end() && !it->second.empty()) {
        args.output_folder = it->second.front();
    }

    bool heat_source_enabled_set = false;
    if (auto it = kv.find("--heat_source_enabled"); it != kv.end() && !it->second.empty()) {
        args.heat_source_enabled = parse_bool_scalar(it->second.front());
        heat_source_enabled_set = true;
    }
    if (auto it = kv.find("--heat_source_min"); it != kv.end()) {
        args.heat_source_min.clear();
        for (const auto& v : it->second) {
            args.heat_source_min.push_back(std::stod(v));
        }
    }
    if (auto it = kv.find("--heat_source_max"); it != kv.end()) {
        args.heat_source_max.clear();
        for (const auto& v : it->second) {
            args.heat_source_max.push_back(std::stod(v));
        }
    }
    if (auto it = kv.find("--heat_source_power_density"); it != kv.end() && !it->second.empty()) {
        args.heat_source_power_density = std::stod(it->second.front());
    }
    if (!heat_source_enabled_set &&
        args.heat_source_min.size() == 3 &&
        args.heat_source_max.size() == 3 &&
        std::abs(args.heat_source_power_density) > 0.0) {
        args.heat_source_enabled = true;
    }

    return args;
}
}  // namespace

// 函数说明：根据文件类型选择 TOML 或旧格式解析路径。
SimulationConfig load_simulation_config(const std::string& path) {
    if (looks_like_toml(path)) {
        return parse_toml_file(path);
    }
    return parse_legacy_file(path);
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
