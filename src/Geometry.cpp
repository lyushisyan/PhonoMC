#include "Geometry.h"

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

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string trim(const std::string& s) {
    const auto b = std::find_if_not(s.begin(), s.end(), [](unsigned char c) { return std::isspace(c); });
    if (b == s.end()) {
        return "";
    }
    const auto e = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) { return std::isspace(c); }).base();
    return std::string(b, e);
}

std::array<long long, 3> quant3(const Vec3& p, double h = 1e-10) {
    return {
        static_cast<long long>(std::llround(p[0] / h)),
        static_cast<long long>(std::llround(p[1] / h)),
        static_cast<long long>(std::llround(p[2] / h))
    };
}

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
T read_binary(std::istream& in) {
    T v {};
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

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

std::vector<Vec3> parse_points(const std::vector<std::string>& tokens, const std::string& name) {
    if (tokens.empty()) {
        return {};
    }
    if ((tokens.size() - 1) % 3 != 0) {
        throw std::runtime_error(name + " points must be declared as type + triplets.");
    }
    std::vector<Vec3> pts;
    pts.reserve((tokens.size() - 1) / 3);
    for (size_t i = 1; i < tokens.size(); i += 3) {
        pts.push_back({std::stod(tokens[i]), std::stod(tokens[i + 1]), std::stod(tokens[i + 2])});
    }
    return pts;
}

std::vector<double> facet_boundary_edge_lengths(const Mesh& mesh, int facet) {
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

Geometry::Geometry(const SimulationConfig& args) {
    build_mesh(args);
    sync_mesh_properties();
    assign_boundary_conditions(args);
    build_periodic_connections(args);
    initialize_subvolumes(args);
    build_subvolume_connections();
    write_geometry_summary(args);

    std::cout << "Geometry initialized: volume=" << volume_
              << ", facets=" << facet_count_
              << ", subvolumes=" << subvolume_count_ << '\n';
}

void Geometry::build_mesh(const SimulationConfig& args) {
    std::vector<Vec3> vertices;
    std::vector<Tri> faces;

    if (args.model == "box") {
        if (args.dimensions.size() != 3) {
            throw std::runtime_error("Box requires 3 dimensions.");
        }
        const Vec3 d {args.dimensions[0], args.dimensions[1], args.dimensions[2]};
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
        if (args.dimensions.size() < 3) {
            throw std::runtime_error("Cylinder requires [L R N].");
        }
        const double L = args.dimensions[0];
        const double R = args.dimensions[1];
        const int N = static_cast<int>(args.dimensions[2]);
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
    }

    mesh_.set_geometry_mesh(std::move(vertices), std::move(faces));
    mesh_.shift_to_origin();
}

void Geometry::sync_mesh_properties() {
    bounds_min_ = mesh_.bounds_min();
    bounds_max_ = mesh_.bounds_max();
    volume_ = mesh_.volume();
    facet_count_ = mesh_.facet_count();
    periodic_pair_.assign(static_cast<size_t>(facet_count_), -1);
    periodic_shift_.assign(static_cast<size_t>(facet_count_), {0.0, 0.0, 0.0});
}

void Geometry::assign_boundary_conditions(const SimulationConfig& args) {
    if (facet_count_ == 0) {
        return;
    }
    const char default_bc = args.boundary_conditions.empty() ? 'P' : args.boundary_conditions.back().empty() ? 'P' : args.boundary_conditions.back()[0];
    facet_boundary_conditions_.assign(facet_count_, default_bc);

    if (!args.boundary_positions.empty()) {
        auto points = parse_points(args.boundary_positions, "boundary_positions");
        if (args.boundary_positions.front() == "relative") {
            const Vec3 ext = sub(bounds_max_, bounds_min_);
            for (auto& p : points) {
                for (int k = 0; k < 3; ++k) {
                    p[k] = bounds_min_[k] + p[k] * ext[k];
                }
            }
        } else if (args.boundary_positions.front() != "absolute") {
            throw std::runtime_error("boundary_positions must start with relative or absolute");
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
            std::fill(roughness_values_.begin(), roughness_values_.end(), args.boundary_values.back());
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
                    roughness_values_[j] = args.boundary_values[val_idx];
                    break;
                }
            }
        }
        ++val_idx;
    }
}

void Geometry::build_periodic_connections(const SimulationConfig& args) {
    connected_facets_.clear();
    if (args.periodic_pair_positions.empty()) {
        return;
    }
    auto points = parse_points(args.periodic_pair_positions, "periodic_pair_positions");
    if (args.periodic_pair_positions.front() == "relative") {
        const Vec3 ext = sub(bounds_max_, bounds_min_);
        for (auto& p : points) {
            for (int k = 0; k < 3; ++k) {
                p[k] = bounds_min_[k] + p[k] * ext[k];
            }
        }
    } else if (args.periodic_pair_positions.front() != "absolute") {
        throw std::runtime_error("periodic_pair_positions must start with relative or absolute");
    }
    if (points.size() % 2 != 0) {
        throw std::runtime_error("periodic_pair_positions must have even number of points.");
    }

    for (size_t i = 0; i < points.size(); i += 2) {
        const int a = mesh_.nearest_facet(points[i]);
        const int b = mesh_.nearest_facet(points[i + 1]);
        connected_facets_.push_back({a, b});
        if (a < 0 || b < 0) {
            throw std::runtime_error("Connected facet mapping failed for periodic_pair_positions pair.");
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

void Geometry::initialize_subvolumes(const SimulationConfig& args) {
    if (args.subvolume_layout.empty()) {
        subvolume_count_ = 1;
        subvolume_centers_ = {mul(add(bounds_min_, bounds_max_), 0.5)};
        subvolume_volumes_ = {volume_};
        return;
    }
    if (args.subvolume_layout.front() == "slice") {
        initialize_slice_subvolumes(args);
    } else if (args.subvolume_layout.front() == "grid") {
        initialize_grid_subvolumes(args);
    } else {
        throw std::runtime_error("Unsupported subvolume layout mode. Use slice or grid.");
    }
}

void Geometry::initialize_slice_subvolumes(const SimulationConfig& args) {
    if (args.subvolume_layout.size() < 3) {
        throw std::runtime_error("slice subvolume layout requires: slice N axis");
    }
    subvolume_count_ = std::stoi(args.subvolume_layout[1]);
    const int axis = std::stoi(args.subvolume_layout[2]);
    if (axis < 0 || axis > 2 || subvolume_count_ <= 0) {
        throw std::runtime_error("Invalid slice settings.");
    }

    subvolume_centers_.assign(static_cast<size_t>(subvolume_count_), mul(add(bounds_min_, bounds_max_), 0.5));
    const double len = bounds_max_[axis] - bounds_min_[axis];
    for (int i = 0; i < subvolume_count_; ++i) {
        const double ratio = (static_cast<double>(i) + 0.5) / static_cast<double>(subvolume_count_);
        subvolume_centers_[i][axis] = bounds_min_[axis] + ratio * len;
    }
    subvolume_volumes_.assign(static_cast<size_t>(subvolume_count_), volume_ / static_cast<double>(subvolume_count_));
}

void Geometry::initialize_grid_subvolumes(const SimulationConfig& args) {
    if (args.subvolume_layout.size() < 4) {
        throw std::runtime_error("grid subvolume layout requires: grid nx ny nz");
    }
    const int nx = std::stoi(args.subvolume_layout[1]);
    const int ny = std::stoi(args.subvolume_layout[2]);
    const int nz = std::stoi(args.subvolume_layout[3]);
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Grid dimensions must be positive.");
    }
    const Vec3 ext = sub(bounds_max_, bounds_min_);
    subvolume_centers_.clear();
    subvolume_centers_.reserve(static_cast<size_t>(nx * ny * nz));
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                const Vec3 p {
                    bounds_min_[0] + (static_cast<double>(ix) + 0.5) * ext[0] / static_cast<double>(nx),
                    bounds_min_[1] + (static_cast<double>(iy) + 0.5) * ext[1] / static_cast<double>(ny),
                    bounds_min_[2] + (static_cast<double>(iz) + 0.5) * ext[2] / static_cast<double>(nz)
                };
                if (mesh_.contains_point(p)) {
                    subvolume_centers_.push_back(p);
                }
            }
        }
    }
    subvolume_count_ = static_cast<int>(subvolume_centers_.size());
    if (subvolume_count_ == 0) {
        throw std::runtime_error("No valid grid subvolume centers found inside mesh.");
    }
    subvolume_volumes_.assign(static_cast<size_t>(subvolume_count_), volume_ / static_cast<double>(subvolume_count_));
}

void Geometry::build_subvolume_connections() {
    subvolume_connections_.clear();
    if (subvolume_count_ <= 1) {
        return;
    }
    for (int i = 0; i < subvolume_count_; ++i) {
        double min_d = std::numeric_limits<double>::max();
        std::vector<std::pair<double, int>> dists;
        dists.reserve(static_cast<size_t>(subvolume_count_ - 1));
        for (int j = 0; j < subvolume_count_; ++j) {
            if (i == j) {
                continue;
            }
            const double d = norm(sub(subvolume_centers_[i], subvolume_centers_[j]));
            dists.push_back({d, j});
            min_d = std::min(min_d, d);
        }
        for (const auto& [d, j] : dists) {
            if (d <= min_d * 1.01) {
                const int a = std::min(i, j);
                const int b = std::max(i, j);
                if (std::find(subvolume_connections_.begin(), subvolume_connections_.end(), std::array<int, 2>{a, b}) == subvolume_connections_.end()) {
                    const Vec3 mid = mul(add(subvolume_centers_[a], subvolume_centers_[b]), 0.5);
                    if (mesh_.contains_point(mid)) {
                        subvolume_connections_.push_back({a, b});
                    } else {
                        const Vec3 dir = sub(subvolume_centers_[b], subvolume_centers_[a]);
                        const auto [_, t, __] = mesh_.trace_boundary_intersection(subvolume_centers_[a], dir);
                        if (std::isinf(t) || t > 1.0) {
                            subvolume_connections_.push_back({a, b});
                        }
                    }
                }
            }
        }
    }
}

void Geometry::write_geometry_summary(const SimulationConfig& args) const {
    if (args.results_base_folder.empty()) {
        return;
    }
    namespace fs = std::filesystem;
    fs::create_directories(args.results_base_folder);
    std::ofstream out(fs::path(args.results_base_folder) / "geometry_summary.txt");
    out << "model=" << args.model << '\n';
    out << "volume=" << volume_ << '\n';
    out << "n_faces=" << mesh_.face_count() << '\n';
    out << "n_facets=" << facet_count_ << '\n';
    out << "n_simplices=" << mesh_.simplex_count() << '\n';
    out << "n_subvols=" << subvolume_count_ << '\n';
    out << "n_connections=" << subvolume_connections_.size() << '\n';
    int npairs = 0;
    for (int i = 0; i < facet_count_; ++i) {
        if (periodic_pair_[i] > i) {
            ++npairs;
        }
    }
    out << "n_periodic_pairs=" << npairs << '\n';
    out << "reservoir_facets=" << reservoir_facets_.size() << '\n';
    out << "rough_facets=" << rough_facets_.size() << '\n';

    std::ofstream centers_out(fs::path(args.results_base_folder) / "subvolume_centers.csv");
    centers_out << "index,x,y,z,volume\n";
    const size_t n = subvolume_centers_.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& c = subvolume_centers_[i];
        const double v = (i < subvolume_volumes_.size()) ? subvolume_volumes_[i] : 0.0;
        centers_out << i << "," << c[0] << "," << c[1] << "," << c[2] << "," << v << "\n";
    }
}

char Geometry::facet_boundary_condition(int facet) const {
    if (facet < 0 || facet >= static_cast<int>(facet_boundary_conditions_.size())) {
        return 'P';
    }
    return facet_boundary_conditions_[facet];
}

double Geometry::reservoir_value_for_facet(int facet, double fallback) const {
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

bool Geometry::has_periodic_pair(int facet) const {
    return facet >= 0 && facet < static_cast<int>(periodic_pair_.size()) && periodic_pair_[facet] >= 0;
}

int Geometry::periodic_pair_facet(int facet) const {
    if (!has_periodic_pair(facet)) {
        return -1;
    }
    return periodic_pair_[facet];
}

std::array<double, 3> Geometry::periodic_shift_for_facet(int facet) const {
    if (!has_periodic_pair(facet)) {
        return {0.0, 0.0, 0.0};
    }
    return periodic_shift_[facet];
}

bool Geometry::is_rough_facet(int facet) const {
    return std::find(rough_facets_.begin(), rough_facets_.end(), facet) != rough_facets_.end();
}

double Geometry::roughness_for_facet(int facet, double fallback) const {
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
