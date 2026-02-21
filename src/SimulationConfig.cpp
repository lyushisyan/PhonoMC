#include "SimulationConfig.h"

#include <algorithm>
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
std::string trim(const std::string& s) {
    const auto first = std::find_if_not(s.begin(), s.end(), [](unsigned char ch) { return std::isspace(ch); });
    if (first == s.end()) {
        return "";
    }
    const auto last = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char ch) { return std::isspace(ch); }).base();
    return std::string(first, last);
}

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
    return s;
}

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

bool is_option(const std::string& tk) {
    return tk.rfind("--", 0) == 0;
}

bool looks_like_toml(const std::string& path) {
    const std::string ext = to_lower(std::filesystem::path(path).extension().string());
    return ext == ".toml";
}

bool is_quoted(const std::string& s) {
    if (s.size() < 2) {
        return false;
    }
    return (s.front() == '"' && s.back() == '"') || (s.front() == '\'' && s.back() == '\'');
}

std::string unquote(const std::string& s) {
    const std::string t = trim(s);
    if (is_quoted(t)) {
        return t.substr(1, t.size() - 2);
    }
    return t;
}

std::string normalise_number_token(std::string s) {
    s = trim(unquote(s));
    s.erase(std::remove(s.begin(), s.end(), '_'), s.end());
    return s;
}

double parse_double_scalar(const std::string& s) {
    const std::string t = normalise_number_token(s);
    if (t.empty()) {
        throw std::runtime_error("Expected numeric value, got empty token.");
    }
    return std::stod(t);
}

int parse_int_scalar(const std::string& s) {
    return static_cast<int>(std::llround(parse_double_scalar(s)));
}

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

std::vector<double> parse_number_array(const std::string& value) {
    std::vector<double> out;
    for (const auto& elem : parse_array_elements(value)) {
        out.push_back(parse_double_scalar(elem));
    }
    return out;
}

std::vector<std::string> parse_string_array(const std::string& value) {
    std::vector<std::string> out;
    for (const auto& elem : parse_array_elements(value)) {
        out.push_back(unquote(elem));
    }
    return out;
}

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

std::string format_number(double x) {
    std::ostringstream oss;
    oss << std::setprecision(17) << x;
    return oss.str();
}

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

SimulationConfig parse_toml_file(const std::string& path) {
    SimulationConfig args;
    args.input_directory = std::filesystem::absolute(std::filesystem::path(path)).parent_path().string();

    std::unordered_map<std::string, std::string> kv;
    parse_toml_assignments(path, kv);

    if (auto v = get_first_value(kv, {"geometry.model", "model"}); v.has_value()) {
        args.model = unquote(*v);
    }
    if (auto v = get_first_value(kv, {"geometry.dimensions", "dimensions"}); v.has_value()) {
        args.dimensions = parse_number_array(*v);
    }
    if (auto v = get_first_value(kv, {"simulation.particle_count", "particle_count", "simulation.particles", "particles"}); v.has_value()) {
        args.particle_count = parse_double_scalar(*v);
    }
    if (auto v = get_first_value(kv, {"simulation.time_step", "time_step", "simulation.timestep", "timestep"}); v.has_value()) {
        args.time_step = parse_double_scalar(*v);
    }
    if (auto v = get_first_value(kv, {"simulation.iterations", "iterations"}); v.has_value()) {
        args.iterations = parse_int_scalar(*v);
    }
    if (auto v = get_first_value(kv, {
            "simulation.compute_thermal_conductivity",
            "simulation.calculate_thermal_conductivity",
            "simulation.calc_kappa",
            "compute_thermal_conductivity",
            "calculate_thermal_conductivity",
            "calc_kappa"}); v.has_value()) {
        args.compute_thermal_conductivity = parse_bool_scalar(*v);
    }

    if (auto v = get_first_value(kv, {"simulation.initial_temperature_profile", "initial_temperature_profile", "simulation.temp_dist", "temp_dist"}); v.has_value()) {
        const std::string t = trim(*v);
        if (!t.empty() && t.front() == '[') {
            args.initial_temperature_profile = parse_string_array(t);
        } else {
            args.initial_temperature_profile = {unquote(t)};
        }
    }

    if (auto v = get_first_value(kv, {"simulation.subvolume_layout", "subvolume_layout", "simulation.subvolumes", "subvolumes"}); v.has_value()) {
        args.subvolume_layout = parse_string_array(*v);
    } else {
        const auto mode_v = get_first_value(kv, {"simulation.subvolumes_mode", "subvolumes_mode"});
        if (mode_v.has_value()) {
            const std::string mode = to_lower(unquote(*mode_v));
            if (mode == "slice") {
                const auto n_v = get_first_value(kv, {"simulation.subvolumes_count", "simulation.subvolumes_n", "subvolumes_count", "subvolumes_n"});
                const auto axis_v = get_first_value(kv, {"simulation.subvolumes_axis", "subvolumes_axis"});
                if (!n_v.has_value() || !axis_v.has_value()) {
                    throw std::runtime_error("TOML slice subvolume layout requires subvolumes_count and subvolumes_axis.");
                }
                args.subvolume_layout = {
                    "slice",
                    std::to_string(parse_int_scalar(*n_v)),
                    std::to_string(parse_int_scalar(*axis_v))
                };
            } else if (mode == "grid") {
                std::vector<int> g(3, 0);
                if (const auto gv = get_first_value(kv, {"simulation.subvolumes_grid", "subvolumes_grid"}); gv.has_value()) {
                    const auto gd = parse_number_array(*gv);
                    if (gd.size() != 3) {
                        throw std::runtime_error("subvolumes_grid must contain 3 values [nx, ny, nz].");
                    }
                    for (int i = 0; i < 3; ++i) {
                        g[static_cast<size_t>(i)] = static_cast<int>(std::llround(gd[static_cast<size_t>(i)]));
                    }
                } else {
                    const auto nx = get_first_value(kv, {"simulation.subvolumes_nx", "subvolumes_nx"});
                    const auto ny = get_first_value(kv, {"simulation.subvolumes_ny", "subvolumes_ny"});
                    const auto nz = get_first_value(kv, {"simulation.subvolumes_nz", "subvolumes_nz"});
                    if (!nx.has_value() || !ny.has_value() || !nz.has_value()) {
                        throw std::runtime_error("TOML grid subvolume layout requires subvolumes_grid or subvolumes_nx/ny/nz.");
                    }
                    g = {parse_int_scalar(*nx), parse_int_scalar(*ny), parse_int_scalar(*nz)};
                }
                args.subvolume_layout = {
                    "grid",
                    std::to_string(g[0]),
                    std::to_string(g[1]),
                    std::to_string(g[2])
                };
            } else {
                throw std::runtime_error("Unsupported subvolumes_mode in TOML: " + mode);
            }
        }
    }

    if (auto v = get_first_value(kv, {"boundary.boundary_conditions", "boundary_conditions", "boundary.bound_cond", "bound_cond"}); v.has_value()) {
        args.boundary_conditions = parse_string_array(*v);
    }
    if (auto v = get_first_value(kv, {"boundary.boundary_values", "boundary_values", "boundary.bound_values", "bound_values"}); v.has_value()) {
        args.boundary_values = parse_number_array(*v);
    }

    if (auto v = get_first_value(kv, {"boundary.boundary_positions", "boundary_positions", "boundary.bound_pos", "bound_pos"}); v.has_value()) {
        args.boundary_positions = parse_string_array(*v);
    } else {
        const auto mode = get_first_value(kv, {"boundary.bound_pos_mode", "bound_pos_mode"});
        const auto pts = get_first_value(kv, {"boundary.bound_pos_points", "bound_pos_points"});
        if (mode.has_value() && pts.has_value()) {
            args.boundary_positions = flatten_points(unquote(*mode), parse_point_array(*pts));
        }
    }

    if (auto v = get_first_value(kv, {"boundary.periodic_pair_positions", "periodic_pair_positions", "boundary.connect_pos", "connect_pos"}); v.has_value()) {
        args.periodic_pair_positions = parse_string_array(*v);
    } else {
        const auto mode = get_first_value(kv, {"boundary.connect_pos_mode", "connect_pos_mode"});
        const auto pts = get_first_value(kv, {"boundary.connect_pos_points", "connect_pos_points"});
        if (mode.has_value() && pts.has_value()) {
            args.periodic_pair_positions = flatten_points(unquote(*mode), parse_point_array(*pts));
        }
    }

    if (auto v = get_first_value(kv, {"io.material_folder", "material_folder", "io.mat_folder", "mat_folder"}); v.has_value()) {
        args.material_folder = unquote(*v);
    }
    if (auto v = get_first_value(kv, {"io.results_base_folder", "results_base_folder", "io.results_folder", "results_folder"}); v.has_value()) {
        args.results_base_folder = unquote(*v);
    }

    bool heat_source_enabled_set = false;
    if (auto v = get_first_value(kv, {"heat_source.enabled", "heat_source_enabled"}); v.has_value()) {
        args.heat_source_enabled = parse_bool_scalar(*v);
        heat_source_enabled_set = true;
    }
    if (auto v = get_first_value(kv, {"heat_source.mode", "heat_source_mode"}); v.has_value()) {
        args.heat_source_mode = to_lower(unquote(*v));
    }
    if (auto v = get_first_value(kv, {"heat_source.min", "heat_source_min"}); v.has_value()) {
        args.heat_source_min = parse_number_array(*v);
    }
    if (auto v = get_first_value(kv, {"heat_source.max", "heat_source_max"}); v.has_value()) {
        args.heat_source_max = parse_number_array(*v);
    }
    if (auto v = get_first_value(kv, {"heat_source.power_density", "heat_source.power_density_wm3", "heat_source_power_density", "heat_source_power_density_wm3"}); v.has_value()) {
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
    if (auto it = kv.find("--dimensions"); it != kv.end() && !it->second.empty()) {
        args.dimensions.clear();
        for (const auto& v : it->second) {
            args.dimensions.push_back(std::stod(v));
        }
    }
    if (auto it = kv.find("--particle_count"); it != kv.end() && !it->second.empty()) {
        args.particle_count = std::stod(it->second.front());
    } else if (auto it2 = kv.find("--particles"); it2 != kv.end() && !it2->second.empty()) {
        args.particle_count = std::stod(it2->second.front());
    }
    if (auto it = kv.find("--time_step"); it != kv.end() && !it->second.empty()) {
        args.time_step = std::stod(it->second.front());
    } else if (auto it2 = kv.find("--timestep"); it2 != kv.end() && !it2->second.empty()) {
        args.time_step = std::stod(it2->second.front());
    }
    if (auto it = kv.find("--iterations"); it != kv.end() && !it->second.empty()) {
        args.iterations = std::stoi(it->second.front());
    }
    if (auto it = kv.find("--compute_thermal_conductivity"); it != kv.end() && !it->second.empty()) {
        args.compute_thermal_conductivity = parse_bool_scalar(it->second.front());
    } else if (auto it2 = kv.find("--calculate_thermal_conductivity"); it2 != kv.end() && !it2->second.empty()) {
        args.compute_thermal_conductivity = parse_bool_scalar(it2->second.front());
    } else if (auto it3 = kv.find("--calc_kappa"); it3 != kv.end() && !it3->second.empty()) {
        args.compute_thermal_conductivity = parse_bool_scalar(it3->second.front());
    }
    if (auto it = kv.find("--subvolume_layout"); it != kv.end()) {
        args.subvolume_layout = it->second;
    } else if (auto it2 = kv.find("--subvolumes"); it2 != kv.end()) {
        args.subvolume_layout = it2->second;
    }
    if (auto it = kv.find("--initial_temperature_profile"); it != kv.end()) {
        args.initial_temperature_profile = it->second;
    } else if (auto it2 = kv.find("--temp_dist"); it2 != kv.end()) {
        args.initial_temperature_profile = it2->second;
    }
    if (auto it = kv.find("--boundary_conditions"); it != kv.end()) {
        args.boundary_conditions = it->second;
    } else if (auto it2 = kv.find("--bound_cond"); it2 != kv.end()) {
        args.boundary_conditions = it2->second;
    }
    if (auto it = kv.find("--boundary_positions"); it != kv.end()) {
        args.boundary_positions = it->second;
    } else if (auto it2 = kv.find("--bound_pos"); it2 != kv.end()) {
        args.boundary_positions = it2->second;
    }
    if (auto it = kv.find("--boundary_values"); it != kv.end()) {
        args.boundary_values.clear();
        for (const auto& v : it->second) {
            args.boundary_values.push_back(std::stod(v));
        }
    } else if (auto it2 = kv.find("--bound_values"); it2 != kv.end()) {
        args.boundary_values.clear();
        for (const auto& v : it2->second) {
            args.boundary_values.push_back(std::stod(v));
        }
    }
    if (auto it = kv.find("--periodic_pair_positions"); it != kv.end()) {
        args.periodic_pair_positions = it->second;
    } else if (auto it2 = kv.find("--connect_pos"); it2 != kv.end()) {
        args.periodic_pair_positions = it2->second;
    }
    if (auto it = kv.find("--material_folder"); it != kv.end() && !it->second.empty()) {
        args.material_folder = it->second.front();
    } else if (auto it2 = kv.find("--mat_folder"); it2 != kv.end() && !it2->second.empty()) {
        args.material_folder = it2->second.front();
    }
    if (auto it = kv.find("--results_base_folder"); it != kv.end() && !it->second.empty()) {
        args.results_base_folder = it->second.front();
    } else if (auto it2 = kv.find("--results_folder"); it2 != kv.end() && !it2->second.empty()) {
        args.results_base_folder = it2->second.front();
    }

    bool heat_source_enabled_set = false;
    if (auto it = kv.find("--heat_source_enabled"); it != kv.end() && !it->second.empty()) {
        args.heat_source_enabled = parse_bool_scalar(it->second.front());
        heat_source_enabled_set = true;
    }
    if (auto it = kv.find("--heat_source_mode"); it != kv.end() && !it->second.empty()) {
        args.heat_source_mode = to_lower(it->second.front());
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
    } else if (auto it2 = kv.find("--heat_source_power_density_wm3"); it2 != kv.end() && !it2->second.empty()) {
        args.heat_source_power_density = std::stod(it2->second.front());
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

SimulationConfig load_simulation_config(const std::string& path) {
    if (looks_like_toml(path)) {
        return parse_toml_file(path);
    }
    return parse_legacy_file(path);
}

std::string create_indexed_results_folder(const std::string& results_base_folder) {
    if (results_base_folder.empty()) {
        throw std::runtime_error("results_base_folder cannot be empty.");
    }

    namespace fs = std::filesystem;
    fs::path base_abs = fs::absolute(fs::path(results_base_folder)).lexically_normal();
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
