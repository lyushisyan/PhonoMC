#include "SimulationDomain.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace {
using Vec3 = std::array<double, 3>;
using Tri = std::array<int, 3>;

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 mul(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// 函数说明：将关键字统一为小写，避免输入大小写差异影响流程。
std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

// 函数说明：清理输入文本两端空白，保证配置与数据解析的稳健性。
std::string trim(const std::string& s) {
    const auto b = std::find_if_not(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); });
    if (b == s.end()) {
        return "";
    }
    const auto e = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) { return std::isspace(c); }).base();
    return std::string(b, e);
}

// 函数说明：将几何坐标量化到容差网格，用于顶点去重键生成。
std::array<long long, 3> quant3(const Vec3& p, double h = 1e-10) {
    return {
        static_cast<long long>(std::llround(p[0] / h)),
        static_cast<long long>(std::llround(p[1] / h)),
        static_cast<long long>(std::llround(p[2] / h))
    };
}

// 函数说明：对顶点做量化去重并返回去重后的顶点索引。
int add_vertex_dedup(
    const Vec3& p,
    std::vector<Vec3>& vertices,
    std::unordered_map<std::string, int>& vmap,
    double h = 1e-10) {
    const auto q = quant3(p, h);
    const std::string k = std::to_string(q[0]) + ":" + std::to_string(q[1]) + ":" + std::to_string(q[2]);
    auto it = vmap.find(k);
    if (it != vmap.end()) {
        return it->second;
    }
    const int idx = static_cast<int>(vertices.size());
    vertices.push_back(p);
    vmap.emplace(k, idx);
    return idx;
}

// 函数说明：读取 OBJ 三角网格并转换为统一顶点/面数据结构。
std::pair<std::vector<Vec3>, std::vector<Tri>> load_obj_mesh(const std::filesystem::path& p) {
    std::ifstream in(p);
    if (!in) {
        throw std::runtime_error("Failed to open OBJ file: " + p.string());
    }
    std::vector<Vec3> vertices;
    std::vector<Tri> faces;
    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') {
            continue;
        }
        std::istringstream iss(t);
        std::string tag;
        iss >> tag;
        if (tag == "v") {
            double x = 0.0, y = 0.0, z = 0.0;
            iss >> x >> y >> z;
            vertices.push_back({x, y, z});
        } else if (tag == "f") {
            std::vector<int> poly;
            std::string tok;
            while (iss >> tok) {
                const auto slash = tok.find('/');
                const std::string idx_s = (slash == std::string::npos) ? tok : tok.substr(0, slash);
                if (idx_s.empty()) {
                    continue;
                }
                int idx = std::stoi(idx_s);
                if (idx < 0) {
                    idx = static_cast<int>(vertices.size()) + idx;
                } else {
                    idx -= 1;
                }
                if (idx < 0 || idx >= static_cast<int>(vertices.size())) {
                    throw std::runtime_error("OBJ face index out of range.");
                }
                poly.push_back(idx);
            }
            if (poly.size() < 3) {
                continue;
            }
            for (size_t i = 1; i + 1 < poly.size(); ++i) {
                faces.push_back({poly[0], poly[i], poly[i + 1]});
            }
        }
    }
    if (vertices.empty() || faces.empty()) {
        throw std::runtime_error("OBJ mesh has no vertices/faces.");
    }
    return {vertices, faces};
}

// 函数说明：读取 ASCII STL 并构建去重后的三角网格。
std::pair<std::vector<Vec3>, std::vector<Tri>> load_ascii_stl_mesh(const std::filesystem::path& p) {
    std::ifstream in(p);
    if (!in) {
        throw std::runtime_error("Failed to open STL file: " + p.string());
    }
    std::vector<Vec3> vertices;
    std::vector<Tri> faces;
    std::unordered_map<std::string, int> vmap;
    std::vector<int> tri_idx;
    tri_idx.reserve(3);

    std::string line;
    while (std::getline(in, line)) {
        const std::string t = trim(line);
        if (t.size() < 6) {
            continue;
        }
        std::istringstream iss(t);
        std::string tag;
        iss >> tag;
        if (tag != "vertex") {
            continue;
        }
        double x = 0.0, y = 0.0, z = 0.0;
        iss >> x >> y >> z;
        tri_idx.push_back(add_vertex_dedup({x, y, z}, vertices, vmap));
        if (tri_idx.size() == 3) {
            faces.push_back({tri_idx[0], tri_idx[1], tri_idx[2]});
            tri_idx.clear();
        }
    }
    if (vertices.empty() || faces.empty()) {
        throw std::runtime_error("ASCII STL mesh has no vertices/faces.");
    }
    return {vertices, faces};
}

template <typename T>
// 函数说明：按二进制类型读取 STL 字段，支撑二进制网格解析。
T read_binary(std::istream& in) {
    T v {};
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

// 函数说明：读取 Binary STL 并转换为内部网格表示。
std::pair<std::vector<Vec3>, std::vector<Tri>> load_binary_stl_mesh(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Failed to open STL file: " + p.string());
    }
    char header[80];
    in.read(header, 80);
    const uint32_t ntri = read_binary<uint32_t>(in);
    std::vector<Vec3> vertices;
    std::vector<Tri> faces;
    std::unordered_map<std::string, int> vmap;
    vertices.reserve(static_cast<size_t>(ntri) * 3);
    faces.reserve(ntri);

    for (uint32_t i = 0; i < ntri; ++i) {
        (void) read_binary<float>(in);
        (void) read_binary<float>(in);
        (void) read_binary<float>(in);
        Vec3 p0 {read_binary<float>(in), read_binary<float>(in), read_binary<float>(in)};
        Vec3 p1 {read_binary<float>(in), read_binary<float>(in), read_binary<float>(in)};
        Vec3 p2 {read_binary<float>(in), read_binary<float>(in), read_binary<float>(in)};
        (void) read_binary<uint16_t>(in);
        const int i0 = add_vertex_dedup(p0, vertices, vmap);
        const int i1 = add_vertex_dedup(p1, vertices, vmap);
        const int i2 = add_vertex_dedup(p2, vertices, vmap);
        faces.push_back({i0, i1, i2});
    }
    if (vertices.empty() || faces.empty()) {
        throw std::runtime_error("Binary STL mesh has no vertices/faces.");
    }
    return {vertices, faces};
}

// 函数说明：依据文件尺寸规则判定 STL 是否为二进制格式。
bool is_probably_binary_stl(const std::filesystem::path& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) {
        return false;
    }
    std::string first5(5, '\0');
    in.read(first5.data(), 5);
    if (!in) {
        return false;
    }
    // Binary STL may also start with "solid"; validate by expected size rule.
    in.seekg(0, std::ios::end);
    const auto sz = static_cast<std::uintmax_t>(in.tellg());
    if (sz < 84) {
        return false;
    }
    in.seekg(80, std::ios::beg);
    const uint32_t ntri = read_binary<uint32_t>(in);
    const std::uintmax_t expected = 84u + static_cast<std::uintmax_t>(ntri) * 50u;
    if (expected == sz) {
        return true;
    }
    return to_lower(first5) != "solid";
}

// 函数说明：按扩展名调度网格加载器，统一支持 OBJ/STL。
std::pair<std::vector<Vec3>, std::vector<Tri>> load_mesh_file(const std::filesystem::path& p) {
    const auto ext = to_lower(p.extension().string());
    if (ext == ".obj") {
        return load_obj_mesh(p);
    }
    if (ext == ".stl") {
        if (is_probably_binary_stl(p)) {
            return load_binary_stl_mesh(p);
        }
        return load_ascii_stl_mesh(p);
    }
    throw std::runtime_error("Unsupported mesh extension: " + ext + " (expected .stl or .obj)");
}

// 函数说明：解析边界/周期点参数并转换为三维坐标列表。
std::vector<Vec3> parse_points(const std::vector<std::string>& tokens, const std::string& name) {
    if (tokens.empty()) {
        return {};
    }
    size_t start = 0;
    try {
        (void) std::stod(tokens.front());
    } catch (...) {
        start = 1;
    }
    if ((tokens.size() - start) % 3 != 0) {
        throw std::runtime_error(name + " points must be declared as triplets.");
    }
    std::vector<Vec3> pts;
    pts.reserve((tokens.size() - start) / 3);
    for (size_t i = start; i < tokens.size(); i += 3) {
        pts.push_back({std::stod(tokens[i]), std::stod(tokens[i + 1]), std::stod(tokens[i + 2])});
    }
    return pts;
}

// 函数说明：提取指定 facet 的边界边长特征，用于边界匹配与诊断。
std::vector<double> facet_boundary_edge_lengths(const SurfaceMesh& mesh, int facet) {
    std::vector<double> out;
    if (facet < 0 || facet >= static_cast<int>(mesh.facet_boundary_edges().size())) {
        return out;
    }
    const auto& eb = mesh.facet_boundary_edges()[facet];
    out.reserve(eb.size());
    for (const int ei : eb) {
        if (ei < 0 || ei >= static_cast<int>(mesh.edges().size())) {
            continue;
        }
        const auto& e = mesh.edges()[ei];
        const auto& v0 = mesh.vertices()[e[0]];
        const auto& v1 = mesh.vertices()[e[1]];
        out.push_back(norm(sub(v0, v1)));
    }
    std::sort(out.begin(), out.end());
    return out;
}
}  // namespace

// 函数说明：构建仿真域：网格、边界、周期配对、控制体与摘要输出。
SimulationDomain::SimulationDomain(const SimulationConfig& args) {
    build_surface_mesh(args);
    sync_surface_mesh_properties();
    assign_boundary_conditions(args);
    build_periodic_connections(args);
    initialize_grids(args);
    build_grid_connections();
    write_domain_summary(args);

    std::cout << "SimulationDomain initialized: volume=" << volume_
              << ", facets=" << facet_count_
              << ", grids=" << grid_count_ << '\n';
}

// 函数说明：根据 box/cylinder/文件模型构建并单位化表面网格。
void SimulationDomain::build_surface_mesh(const SimulationConfig& args) {
    std::vector<Vec3> vertices;
    std::vector<Tri> faces;

    if (args.model == "box") {
        if (args.sizes.size() != 3) {
            throw std::runtime_error("Box requires 3 sizes.");
        }
        const Vec3 d {args.sizes[0], args.sizes[1], args.sizes[2]};
        vertices = {
            Vec3{0, 0, 0}, Vec3{0, 0, 1}, Vec3{0, 1, 1}, Vec3{0, 1, 0},
            Vec3{1, 0, 0}, Vec3{1, 0, 1}, Vec3{1, 1, 1}, Vec3{1, 1, 0}
        };
        for (auto& v : vertices) {
            v[0] *= d[0];
            v[1] *= d[1];
            v[2] *= d[2];
        }
        faces = {
            Tri{0, 1, 2}, Tri{0, 2, 3}, Tri{4, 5, 6}, Tri{4, 6, 7},
            Tri{0, 4, 5}, Tri{0, 5, 1}, Tri{3, 7, 6}, Tri{3, 6, 2},
            Tri{0, 4, 7}, Tri{0, 7, 3}, Tri{1, 5, 6}, Tri{1, 6, 2}
        };
    } else if (args.model == "cylinder") {
        if (args.sizes.size() < 3) {
            throw std::runtime_error("Cylinder requires [L R N].");
        }
        const double L = args.sizes[0];
        const double R = args.sizes[1];
        const int N = static_cast<int>(args.sizes[2]);
        if (N < 3) {
            throw std::runtime_error("Cylinder segment count N must be >= 3.");
        }
        vertices.push_back({0.0, 0.0, 0.0});
        for (int i = 0; i < N; ++i) {
            const double a = (2.0 * M_PI * i) / static_cast<double>(N);
            vertices.push_back({R * std::cos(a), R * std::sin(a), 0.0});
        }
        vertices.push_back({0.0, 0.0, L});
        for (int i = 0; i < N; ++i) {
            const double a = (2.0 * M_PI * i) / static_cast<double>(N);
            vertices.push_back({R * std::cos(a), R * std::sin(a), L});
        }
        for (int i = 1; i <= N; ++i) {
            const int j = (i == N) ? 1 : (i + 1);
            faces.push_back({0, i, j});
        }
        for (int i = 1; i <= N; ++i) {
            const int j = (i == N) ? 1 : (i + 1);
            const int io = i + (N + 1);
            const int jo = j + (N + 1);
            faces.push_back({i, io, jo});
            faces.push_back({i, j, jo});
        }
        const int top_center = N + 1;
        for (int i = 1; i <= N; ++i) {
            const int j = (i == N) ? 1 : (i + 1);
            faces.push_back({top_center, j + (N + 1), i + (N + 1)});
        }
    } else {
        std::filesystem::path mesh_path(args.model);
        if (!mesh_path.is_absolute()) {
            std::filesystem::path p1 = std::filesystem::current_path() / mesh_path;
            std::filesystem::path p2 = std::filesystem::path(args.input_directory) / mesh_path;
            if (std::filesystem::exists(p1)) {
                mesh_path = p1;
            } else {
                mesh_path = p2;
            }
        }
        mesh_path = mesh_path.lexically_normal();
        if (!std::filesystem::exists(mesh_path)) {
            throw std::runtime_error("Model file does not exist: " + mesh_path.string());
        }
        auto loaded = load_mesh_file(mesh_path);
        vertices = std::move(loaded.first);
        faces = std::move(loaded.second);
        // STL geometry is provided in nm; convert to internal Angstrom.
        if (to_lower(mesh_path.extension().string()) == ".stl") {
            for (auto& v : vertices) {
                v[0] *= 10.0;
                v[1] *= 10.0;
                v[2] *= 10.0;
            }
        }
    }

    mesh_.set_surface_mesh_data(std::move(vertices), std::move(faces));
    mesh_.shift_to_origin();
}

// 函数说明：同步网格边界框、体积和 facet 数等几何缓存。
void SimulationDomain::sync_surface_mesh_properties() {
    bounds_min_ = mesh_.bounds_min();
    bounds_max_ = mesh_.bounds_max();
    volume_ = mesh_.volume();
    facet_count_ = mesh_.facet_count();
    periodic_pair_.assign(static_cast<size_t>(facet_count_), -1);
    periodic_shift_.assign(static_cast<size_t>(facet_count_), {0.0, 0.0, 0.0});
}

// 函数说明：将输入边界条件映射到 facet，并识别热库/粗糙边界。
void SimulationDomain::assign_boundary_conditions(const SimulationConfig& args) {
    if (facet_count_ == 0) {
        return;
    }
    const char default_bc = args.boundary_conditions.empty() ? 'P' : args.boundary_conditions.back().empty() ? 'P' : args.boundary_conditions.back()[0];
    facet_boundary_conditions_.assign(facet_count_, default_bc);

    if (!args.boundary_position.empty()) {
        auto points = parse_points(args.boundary_position, "boundary_position");
        const Vec3 ext = sub(bounds_max_, bounds_min_);
        for (auto& p : points) {
            for (int k = 0; k < 3; ++k) {
                p[k] = bounds_min_[k] + p[k] * ext[k];
            }
        }

        boundary_facets_.clear();
        for (size_t i = 0; i < points.size(); ++i) {
            const int best_f = mesh_.nearest_facet(points[i]);
            if (best_f < 0) {
                continue;
            }
            boundary_facets_.push_back(best_f);
            if (i < args.boundary_conditions.size() && !args.boundary_conditions[i].empty()) {
                facet_boundary_conditions_[best_f] = args.boundary_conditions[i][0];
            }
        }
    }

    reservoir_facets_.clear();
    rough_facets_.clear();
    for (int f = 0; f < facet_count_; ++f) {
        if (facet_boundary_conditions_[f] == 'T' || facet_boundary_conditions_[f] == 'F') {
            reservoir_facets_.push_back(f);
        } else if (facet_boundary_conditions_[f] == 'R') {
            rough_facets_.push_back(f);
        }
    }

    reservoir_values_.assign(reservoir_facets_.size(), std::numeric_limits<double>::quiet_NaN());
    roughness_values_.assign(rough_facets_.size(), std::numeric_limits<double>::quiet_NaN());

    if (!args.boundary_values.empty()) {
        if (default_bc == 'T' || default_bc == 'F') {
            std::fill(reservoir_values_.begin(), reservoir_values_.end(), args.boundary_values.back());
        } else if (default_bc == 'R') {
            // Roughness is provided in nm in input; convert to Angstrom internally.
            std::fill(roughness_values_.begin(), roughness_values_.end(), args.boundary_values.back() * 10.0);
        }
    }

    int val_idx = 0;
    for (size_t i = 0; i < boundary_facets_.size(); ++i) {
        const int f = boundary_facets_[i];
        const char bc = facet_boundary_conditions_[f];
        if (bc == 'P') {
            continue;
        }
        if (val_idx >= static_cast<int>(args.boundary_values.size())) {
            break;
        }
        if (bc == 'T' || bc == 'F') {
            for (size_t j = 0; j < reservoir_facets_.size(); ++j) {
                if (reservoir_facets_[j] == f) {
                    reservoir_values_[j] = args.boundary_values[val_idx];
                    break;
                }
            }
        } else if (bc == 'R') {
            for (size_t j = 0; j < rough_facets_.size(); ++j) {
                if (rough_facets_[j] == f) {
                    // Roughness is provided in nm in input; convert to Angstrom internally.
                    roughness_values_[j] = args.boundary_values[val_idx] * 10.0;
                    break;
                }
            }
        }
        ++val_idx;
    }
}

// 函数说明：建立周期边界配对关系与平移向量。
void SimulationDomain::build_periodic_connections(const SimulationConfig& args) {
    connected_facets_.clear();
    if (args.periodic_pair.empty()) {
        return;
    }
    auto points = parse_points(args.periodic_pair, "periodic_pair");
    const Vec3 ext = sub(bounds_max_, bounds_min_);
    for (auto& p : points) {
        for (int k = 0; k < 3; ++k) {
            p[k] = bounds_min_[k] + p[k] * ext[k];
        }
    }
    if (points.size() % 2 != 0) {
        throw std::runtime_error("periodic_pair must have even number of points.");
    }

    for (size_t i = 0; i < points.size(); i += 2) {
        const int a = mesh_.nearest_facet(points[i]);
        const int b = mesh_.nearest_facet(points[i + 1]);
        connected_facets_.push_back({a, b});
        if (a < 0 || b < 0) {
            throw std::runtime_error("Connected facet mapping failed for periodic_pair pair.");
        }
        if (facet_boundary_conditions_[a] != 'P' || facet_boundary_conditions_[b] != 'P') {
            throw std::runtime_error("Connected facets must both be periodic ('P').");
        }
        if ((periodic_pair_[a] >= 0 && periodic_pair_[a] != b) ||
            (periodic_pair_[b] >= 0 && periodic_pair_[b] != a)) {
            throw std::runtime_error("Facet already paired with a different periodic connection.");
        }
        const double ndot = dot(mesh_.facet_normals()[a], mesh_.facet_normals()[b]);
        if (ndot > -0.98) {
            throw std::runtime_error("Connected facets normals do not match opposite direction.");
        }
        const double aa = mesh_.facet_areas()[a];
        const double ab = mesh_.facet_areas()[b];
        const double area_rel = std::abs(aa - ab) / std::max({1.0, aa, ab});
        if (area_rel > 1e-3) {
            throw std::runtime_error("Connected facets area mismatch.");
        }
        const auto la = facet_boundary_edge_lengths(mesh_, a);
        const auto lb = facet_boundary_edge_lengths(mesh_, b);
        if (la.size() != lb.size()) {
            throw std::runtime_error("Connected facets boundary edge count mismatch.");
        }
        for (size_t k = 0; k < la.size(); ++k) {
            const double rel = std::abs(la[k] - lb[k]) / std::max({1.0, la[k], lb[k]});
            if (rel > 1e-3) {
                throw std::runtime_error("Connected facets boundary edge-length spectrum mismatch.");
            }
        }

        // Register periodic translation mapping (facet a -> b and b -> a).
        const Vec3 shift_ab = sub(mesh_.facet_centroids()[b], mesh_.facet_centroids()[a]);
        const Vec3 shift_ba = sub(mesh_.facet_centroids()[a], mesh_.facet_centroids()[b]);
        periodic_pair_[a] = b;
        periodic_pair_[b] = a;
        periodic_shift_[a] = shift_ab;
        periodic_shift_[b] = shift_ba;
    }
}

// 函数说明：初始化控制体网格布局入口（当前为规则网格模式）。
void SimulationDomain::initialize_grids(const SimulationConfig& args) {
    if (args.grid_layout.front() == "grid") {
        initialize_grid_cells(args);
    } else {
        throw std::runtime_error("Unsupported grid layout mode. Use grid.");
    }
}

// 函数说明：在几何域内生成有效网格中心与体积权重。
void SimulationDomain::initialize_grid_cells(const SimulationConfig& args) {
    if (args.grid_layout.size() < 4) {
        throw std::runtime_error("grid layout requires: grid nx ny nz");
    }
    const int nx = std::stoi(args.grid_layout[1]);
    const int ny = std::stoi(args.grid_layout[2]);
    const int nz = std::stoi(args.grid_layout[3]);
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Grid sizes must be positive.");
    }
    const Vec3 ext = sub(bounds_max_, bounds_min_);
    grid_centers_.clear();
    grid_centers_.reserve(static_cast<size_t>(nx * ny * nz));
    const bool is_box = (args.model == "box");
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                const Vec3 p {
                    bounds_min_[0] + (static_cast<double>(ix) + 0.5) * ext[0] / static_cast<double>(nx),
                    bounds_min_[1] + (static_cast<double>(iy) + 0.5) * ext[1] / static_cast<double>(ny),
                    bounds_min_[2] + (static_cast<double>(iz) + 0.5) * ext[2] / static_cast<double>(nz)
                };
                if (is_box || mesh_.contains_point(p)) {
                    grid_centers_.push_back(p);
                }
            }
        }
    }
    grid_count_ = static_cast<int>(grid_centers_.size());
    if (grid_count_ == 0) {
        throw std::runtime_error("No valid grid centers found inside mesh.");
    }
    grid_volumes_.assign(static_cast<size_t>(grid_count_), volume_ / static_cast<double>(grid_count_));
}

// 函数说明：基于邻近关系与几何可达性构建网格连接图。
void SimulationDomain::build_grid_connections() {
    grid_connections_.clear();
    if (grid_count_ <= 1) {
        return;
    }
    for (int i = 0; i < grid_count_; ++i) {
        double min_d = std::numeric_limits<double>::max();
        std::vector<std::pair<double, int>> dists;
        dists.reserve(static_cast<size_t>(grid_count_ - 1));
        for (int j = 0; j < grid_count_; ++j) {
            if (i == j) {
                continue;
            }
            const double d = norm(sub(grid_centers_[i], grid_centers_[j]));
            dists.push_back({d, j});
            min_d = std::min(min_d, d);
        }
        for (const auto& [d, j] : dists) {
            if (d <= min_d * 1.01) {
                const int a = std::min(i, j);
                const int b = std::max(i, j);
                if (std::find(grid_connections_.begin(), grid_connections_.end(), std::array<int, 2>{a, b}) == grid_connections_.end()) {
                    const Vec3 mid = mul(add(grid_centers_[a], grid_centers_[b]), 0.5);
                    if (mesh_.contains_point(mid)) {
                        grid_connections_.push_back({a, b});
                    } else {
                        const Vec3 dir = sub(grid_centers_[b], grid_centers_[a]);
                        const auto [_, t, __] = mesh_.trace_boundary_intersection(grid_centers_[a], dir);
                        if (std::isinf(t) || t > 1.0) {
                            grid_connections_.push_back({a, b});
                        }
                    }
                }
            }
        }
    }
}

// 函数说明：输出几何、边界、网格与输入信息摘要文件。
void SimulationDomain::write_domain_summary(const SimulationConfig& args) const {
    if (args.output_folder.empty()) {
        return;
    }
    namespace fs = std::filesystem;
    fs::create_directories(args.output_folder);
    std::ofstream out(fs::path(args.output_folder) / "summary.txt");
    auto join_strings = [](const std::vector<std::string>& v) -> std::string {
        if (v.empty()) {
            return "[]";
        }
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) {
                oss << ", ";
            }
            oss << v[i];
        }
        oss << "]";
        return oss.str();
    };
    auto join_numbers = [](const std::vector<double>& v) -> std::string {
        if (v.empty()) {
            return "[]";
        }
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) {
                oss << ", ";
            }
            oss << v[i];
        }
        oss << "]";
        return oss.str();
    };
    auto point_count = [](const std::vector<std::string>& v) -> int {
        if (v.empty()) {
            return 0;
        }
        size_t start = 0;
        try {
            (void) std::stod(v.front());
        } catch (...) {
            start = 1;
        }
        if ((v.size() - start) % 3 != 0) {
            return 0;
        }
        return static_cast<int>((v.size() - start) / 3);
    };

    int npairs = 0;
    for (int i = 0; i < facet_count_; ++i) {
        if (periodic_pair_[i] > i) {
            ++npairs;
        }
    }

    out << "[input]\n";
    out << "model = " << args.model << '\n';
    out << "sizes_nm = [" << (args.sizes.size() > 0 ? args.sizes[0] / 10.0 : 0.0)
        << ", " << (args.sizes.size() > 1 ? args.sizes[1] / 10.0 : 0.0)
        << ", " << (args.sizes.size() > 2 ? args.sizes[2] / 10.0 : 0.0) << "]\n";
    out << "sizes_angstrom = [" << (args.sizes.size() > 0 ? args.sizes[0] : 0.0)
        << ", " << (args.sizes.size() > 1 ? args.sizes[1] : 0.0)
        << ", " << (args.sizes.size() > 2 ? args.sizes[2] : 0.0) << "]\n";
    out << "particle_count = " << args.particle_count << '\n';
    out << "time_step = " << args.time_step << '\n';
    out << "iterations = " << args.iterations << '\n';
    out << "compute_kappa = " << (args.compute_kappa ? "true" : "false") << '\n';
    out << "profile_timers = " << (args.profile_timers ? "true" : "false") << '\n';
    out << "initial_temperature = " << join_strings(args.initial_temperature) << '\n';
    if (args.grid_layout.size() >= 4) {
        out << "grid_xyz = [" << args.grid_layout[1] << ", " << args.grid_layout[2] << ", " << args.grid_layout[3] << "]\n";
    } else {
        out << "grid_xyz = []\n";
    }
    out << "boundary_conditions = " << join_strings(args.boundary_conditions) << '\n';
    out << "boundary_values = " << join_numbers(args.boundary_values) << '\n';
    out << "boundary_position_point_count = " << point_count(args.boundary_position) << '\n';
    out << "boundary_position_tokens = " << join_strings(args.boundary_position) << '\n';
    out << "periodic_pair_point_count = " << point_count(args.periodic_pair) << '\n';
    out << "periodic_pair_tokens = " << join_strings(args.periodic_pair) << '\n';
    out << "material_folder = " << args.material_folder << '\n';
    out << "output_folder = " << args.output_folder << '\n';
    out << "heat_source_enabled = " << (args.heat_source_enabled ? "true" : "false") << '\n';
    out << "heat_source_min = " << join_numbers(args.heat_source_min) << '\n';
    out << "heat_source_max = " << join_numbers(args.heat_source_max) << '\n';
    out << "heat_source_power_density = " << args.heat_source_power_density << '\n';
    out << '\n';

    out << "[geometry]\n";
    out << "volume = " << volume_ << '\n';
    out << "bounds_min = [" << bounds_min_[0] << ", " << bounds_min_[1] << ", " << bounds_min_[2] << "]\n";
    out << "bounds_max = [" << bounds_max_[0] << ", " << bounds_max_[1] << ", " << bounds_max_[2] << "]\n";
    out << "n_faces = " << mesh_.face_count() << '\n';
    out << "n_facets = " << facet_count_ << '\n';
    out << "n_simplices = " << mesh_.simplex_count() << '\n';
    out << '\n';

    out << "[grid]\n";
    out << "count = " << grid_count_ << '\n';
    out << "connections = " << grid_connections_.size() << '\n';
    out << '\n';

    out << "[boundary]\n";
    out << "reservoir_facets = " << reservoir_facets_.size() << '\n';
    out << "rough_facets = " << rough_facets_.size() << '\n';
    out << "periodic_pairs = " << npairs << '\n';
    out << '\n';

    out << "[files]\n";
    out << "grid_centers_csv = " << (fs::path(args.output_folder) / "grid_centers.csv").string() << '\n';

    std::ofstream centers_out(fs::path(args.output_folder) / "grid_centers.csv");
    centers_out << "index,x_nm,y_nm,z_nm,volume_nm3\n";
    const size_t n = grid_centers_.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& c = grid_centers_[i];
        const double v = (i < grid_volumes_.size()) ? grid_volumes_[i] : 0.0;
        // Internal length unit is Angstrom; export centers in nm for user-facing outputs.
        centers_out << i << ","
                    << (c[0] / 10.0) << ","
                    << (c[1] / 10.0) << ","
                    << (c[2] / 10.0) << ","
                    << (v / 1000.0) << "\n";
    }
}

// 函数说明：查询指定 facet 的边界类型（T/P/R/F 等）。
char SimulationDomain::facet_boundary_condition(int facet) const {
    if (facet < 0 || facet >= static_cast<int>(facet_boundary_conditions_.size())) {
        return 'P';
    }
    return facet_boundary_conditions_[facet];
}

// 函数说明：查询热库 facet 的边界温度值。
double SimulationDomain::reservoir_value_for_facet(int facet, double fallback) const {
    for (size_t i = 0; i < reservoir_facets_.size(); ++i) {
        if (reservoir_facets_[i] == facet) {
            if (i < reservoir_values_.size() && !std::isnan(reservoir_values_[i])) {
                return reservoir_values_[i];
            }
            return fallback;
        }
    }
    return fallback;
}

// 函数说明：判断 facet 是否存在周期配对面。
bool SimulationDomain::has_periodic_pair(int facet) const {
    return facet >= 0 && facet < static_cast<int>(periodic_pair_.size()) && periodic_pair_[facet] >= 0;
}

// 函数说明：返回 facet 对应的周期配对面索引。
int SimulationDomain::periodic_pair_facet(int facet) const {
    if (!has_periodic_pair(facet)) {
        return -1;
    }
    return periodic_pair_[facet];
}

// 函数说明：返回周期穿越时粒子位置平移向量。
std::array<double, 3> SimulationDomain::periodic_shift_for_facet(int facet) const {
    if (!has_periodic_pair(facet)) {
        return {0.0, 0.0, 0.0};
    }
    return periodic_shift_[facet];
}

// 函数说明：判断 facet 是否为粗糙散射边界。
bool SimulationDomain::is_rough_facet(int facet) const {
    return std::find(rough_facets_.begin(), rough_facets_.end(), facet) != rough_facets_.end();
}

// 函数说明：查询 facet 粗糙度参数，用于镜面率计算。
double SimulationDomain::roughness_for_facet(int facet, double fallback) const {
    for (size_t i = 0; i < rough_facets_.size(); ++i) {
        if (rough_facets_[i] == facet) {
            if (i < roughness_values_.size() && !std::isnan(roughness_values_[i])) {
                return roughness_values_[i];
            }
            return fallback;
        }
    }
    return fallback;
}
