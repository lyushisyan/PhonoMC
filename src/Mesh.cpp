#include "Mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {
using Vec3 = Mesh::Vec3;

Vec3 add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 mul(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 normalize(const Vec3& a) {
    const double n = norm(a);
    if (n <= 0.0) {
        return {0.0, 0.0, 0.0};
    }
    return mul(a, 1.0 / n);
}
Vec3 min_vec(const Vec3& a, const Vec3& b) {
    return {std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2])};
}
Vec3 max_vec(const Vec3& a, const Vec3& b) {
    return {std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2])};
}

std::array<std::array<double, 3>, 3> inverse3x3(const std::array<std::array<double, 3>, 3>& m, bool& ok) {
    const double a = m[0][0], b = m[0][1], c = m[0][2];
    const double d = m[1][0], e = m[1][1], f = m[1][2];
    const double g = m[2][0], h = m[2][1], i = m[2][2];
    const double A = e * i - f * h;
    const double B = -(d * i - f * g);
    const double C = d * h - e * g;
    const double D = -(b * i - c * h);
    const double E = a * i - c * g;
    const double F = -(a * h - b * g);
    const double G = b * f - c * e;
    const double H = -(a * f - c * d);
    const double I = a * e - b * d;
    const double det = a * A + b * B + c * C;
    if (std::abs(det) < 1e-14) {
        ok = false;
        return {};
    }
    const double id = 1.0 / det;
    ok = true;
    return {{
        {{A * id, D * id, G * id}},
        {{B * id, E * id, H * id}},
        {{C * id, F * id, I * id}}
    }};
}

Vec3 matvec3(const std::array<std::array<double, 3>, 3>& m, const Vec3& v) {
    return {
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
    };
}

std::string quant_key(const Vec3& p, double h) {
    auto q = [h](double x) -> long long {
        return static_cast<long long>(std::llround(x / h));
    };
    return std::to_string(q(p[0])) + ":" + std::to_string(q(p[1])) + ":" + std::to_string(q(p[2]));
}

Vec3 closest_point_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 ab = sub(b, a);
    const Vec3 ac = sub(c, a);
    const Vec3 ap = sub(p, a);

    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        return a;
    }

    const Vec3 bp = sub(p, b);
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        return b;
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        return add(a, mul(ab, v));
    }

    const Vec3 cp = sub(p, c);
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        return c;
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        return add(a, mul(ac, w));
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const Vec3 bc = sub(c, b);
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, mul(bc, w));
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    return add(add(a, mul(ab, v)), mul(ac, w));
}

bool ray_intersects_triangle(const Vec3& orig, const Vec3& dir, const Vec3& a, const Vec3& b, const Vec3& c, double& t) {
    constexpr double eps = 1e-10;
    const Vec3 e1 = sub(b, a);
    const Vec3 e2 = sub(c, a);
    const Vec3 pvec = cross(dir, e2);
    const double det = dot(e1, pvec);
    if (std::abs(det) < eps) {
        return false;
    }
    const double inv_det = 1.0 / det;
    const Vec3 tvec = sub(orig, a);
    const double u = dot(tvec, pvec) * inv_det;
    if (u < -eps || u > 1.0 + eps) {
        return false;
    }
    const Vec3 qvec = cross(tvec, e1);
    const double v = dot(dir, qvec) * inv_det;
    if (v < -eps || (u + v) > 1.0 + eps) {
        return false;
    }
    t = dot(e2, qvec) * inv_det;
    return t > eps;
}
}  // namespace

Mesh::Mesh(std::vector<Vec3> vertices, std::vector<Tri> faces) {
    set_geometry_mesh(std::move(vertices), std::move(faces));
}

void Mesh::set_geometry_mesh(std::vector<Vec3> vertices, std::vector<Tri> faces) {
    vertices_ = std::move(vertices);
    faces_ = std::move(faces);
    if (vertices_.empty() || faces_.empty()) {
        throw std::runtime_error("Mesh requires non-empty vertices and faces.");
    }
    clear_tetrahedra_cache();
    rebuild_cached_properties();
}

void Mesh::shift_to_origin() {
    compute_bounding_box();
    for (auto& v : vertices_) {
        v = sub(v, bounds_min_);
    }
    rebuild_cached_properties();
}

void Mesh::rebuild_cached_properties() {
    clear_tetrahedra_cache();
    compute_bounding_box();
    compute_edge_list();
    compute_face_metrics();
    compute_face_neighbors();
    compute_facet_groups();
    compute_facet_neighbors();
    compute_interface_faces();
    compute_enclosed_volume();
}

void Mesh::clear_tetrahedra_cache() {
    simplices_.clear();
    simplices_ready_ = false;
}

void Mesh::compute_bounding_box() {
    bounds_min_ = vertices_.front();
    bounds_max_ = vertices_.front();
    for (const auto& v : vertices_) {
        bounds_min_ = min_vec(bounds_min_, v);
        bounds_max_ = max_vec(bounds_max_, v);
    }
}

void Mesh::compute_edge_list() {
    std::map<std::pair<int, int>, int> edge_map;
    for (const auto& f : faces_) {
        const std::array<std::array<int, 2>, 3> e {
            std::array<int, 2> {std::min(f[0], f[1]), std::max(f[0], f[1])},
            std::array<int, 2> {std::min(f[0], f[2]), std::max(f[0], f[2])},
            std::array<int, 2> {std::min(f[1], f[2]), std::max(f[1], f[2])}
        };
        for (const auto& seg : e) {
            edge_map.emplace(std::make_pair(seg[0], seg[1]), static_cast<int>(edge_map.size()));
        }
    }
    edges_.assign(edge_map.size(), {0, 0});
    for (const auto& [k, idx] : edge_map) {
        edges_[idx] = {k.first, k.second};
    }

    face_edges_.assign(faces_.size(), {0, 0, 0});
    edges_faces_.assign(edges_.size(), {});
    for (size_t i = 0; i < faces_.size(); ++i) {
        const auto& f = faces_[i];
        const std::array<std::pair<int, int>, 3> e {
            std::make_pair(std::min(f[0], f[1]), std::max(f[0], f[1])),
            std::make_pair(std::min(f[0], f[2]), std::max(f[0], f[2])),
            std::make_pair(std::min(f[1], f[2]), std::max(f[1], f[2]))
        };
        for (int j = 0; j < 3; ++j) {
            const int ei = edge_map.at(e[j]);
            face_edges_[i][j] = ei;
            edges_faces_[ei].push_back(static_cast<int>(i));
        }
    }
}

void Mesh::compute_face_metrics() {
    face_normals_.assign(faces_.size(), Vec3 {0.0, 0.0, 0.0});
    face_centroids_.assign(faces_.size(), Vec3 {0.0, 0.0, 0.0});
    face_areas_.assign(faces_.size(), 0.0);
    face_k_.assign(faces_.size(), 0.0);
    face_bounds_.assign(faces_.size(), {Vec3 {0, 0, 0}, Vec3 {0, 0, 0}});

    for (size_t i = 0; i < faces_.size(); ++i) {
        const auto& f = faces_[i];
        const Vec3& a = vertices_[f[0]];
        const Vec3& b = vertices_[f[1]];
        const Vec3& c = vertices_[f[2]];
        const Vec3 n = cross(sub(b, a), sub(c, a));
        const double area2 = norm(n);
        face_normals_[i] = (area2 > 0.0) ? mul(n, 1.0 / area2) : Vec3 {0, 0, 0};
        face_areas_[i] = area2 * 0.5;
        face_centroids_[i] = {(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0, (a[2] + b[2] + c[2]) / 3.0};
        face_k_[i] = -dot(face_normals_[i], a);
        face_bounds_[i] = {min_vec(min_vec(a, b), c), max_vec(max_vec(a, b), c)};
    }
}

void Mesh::compute_face_neighbors() {
    face_adjacency_.clear();
    for (const auto& ef : edges_faces_) {
        if (ef.size() < 2) {
            continue;
        }
        for (size_t i = 0; i < ef.size(); ++i) {
            for (size_t j = i + 1; j < ef.size(); ++j) {
                const int a = std::min(ef[i], ef[j]);
                const int b = std::max(ef[i], ef[j]);
                face_adjacency_.push_back({a, b});
            }
        }
    }
    std::sort(face_adjacency_.begin(), face_adjacency_.end(), [](const auto& a, const auto& b) {
        if (a[0] != b[0]) {
            return a[0] < b[0];
        }
        return a[1] < b[1];
    });
    face_adjacency_.erase(std::unique(face_adjacency_.begin(), face_adjacency_.end()), face_adjacency_.end());
}

void Mesh::compute_facet_groups() {
    face_to_facet_.assign(faces_.size(), -1);
    facets_.clear();

    struct FacetAccum {
        Vec3 normal {};
        Vec3 centroid_weighted {};
        double area = 0.0;
        std::vector<int> face_indices;
    };
    std::unordered_map<std::string, FacetAccum> acc;

    auto key_of = [](const Vec3& n, double d) {
        auto q = [](double x) { return static_cast<long long>(std::llround(x * 1e6)); };
        return std::to_string(q(n[0])) + ":" + std::to_string(q(n[1])) + ":" + std::to_string(q(n[2])) + ":" + std::to_string(q(d));
    };

    const Vec3 center = mul(add(bounds_min_, bounds_max_), 0.5);
    for (size_t i = 0; i < faces_.size(); ++i) {
        Vec3 n = face_normals_[i];
        if (dot(n, sub(face_centroids_[i], center)) < 0.0) {
            n = mul(n, -1.0);
        }
        const double d = dot(n, vertices_[faces_[i][0]]);
        auto& a = acc[key_of(n, d)];
        a.normal = add(a.normal, n);
        a.centroid_weighted = add(a.centroid_weighted, mul(face_centroids_[i], face_areas_[i]));
        a.area += face_areas_[i];
        a.face_indices.push_back(static_cast<int>(i));
    }

    facets_.reserve(acc.size());
    int idx = 0;
    for (auto& [_, a] : acc) {
        Facet f;
        f.normal = normalize(a.normal);
        f.area = a.area;
        f.centroid = (a.area > 0.0) ? mul(a.centroid_weighted, 1.0 / a.area) : a.centroid_weighted;
        f.faces = std::move(a.face_indices);
        for (const int fi : f.faces) {
            face_to_facet_[fi] = idx;
        }
        facets_.push_back(std::move(f));
        ++idx;
    }

    facets_edges_.clear();
    facets_boundary_edges_.clear();
    facet_areas_.clear();
    facet_normals_.clear();
    facet_centroids_.clear();
    facets_edges_.reserve(facets_.size());
    facets_boundary_edges_.reserve(facets_.size());
    facet_areas_.reserve(facets_.size());
    facet_normals_.reserve(facets_.size());
    facet_centroids_.reserve(facets_.size());
    for (const auto& f : facets_) {
        std::unordered_map<int, int> edge_count;
        for (const int fi : f.faces) {
            for (const int ei : face_edges_[fi]) {
                edge_count[ei] += 1;
            }
        }
        std::vector<int> all_edges;
        std::vector<int> boundary_edges;
        all_edges.reserve(edge_count.size());
        boundary_edges.reserve(edge_count.size());
        for (const auto& [ei, cnt] : edge_count) {
            all_edges.push_back(ei);
            if (cnt == 1) {
                boundary_edges.push_back(ei);
            }
        }
        std::sort(all_edges.begin(), all_edges.end());
        std::sort(boundary_edges.begin(), boundary_edges.end());
        facets_edges_.push_back(std::move(all_edges));
        facets_boundary_edges_.push_back(std::move(boundary_edges));
        facet_areas_.push_back(f.area);
        facet_normals_.push_back(f.normal);
        facet_centroids_.push_back(f.centroid);
    }
}

void Mesh::compute_facet_neighbors() {
    facets_adjacency_.clear();
    for (const auto& adj : face_adjacency_) {
        const int fa = face_to_facet_[adj[0]];
        const int fb = face_to_facet_[adj[1]];
        if (fa < 0 || fb < 0 || fa == fb) {
            continue;
        }
        const int a = std::min(fa, fb);
        const int b = std::max(fa, fb);
        facets_adjacency_.push_back({a, b});
    }
    std::sort(facets_adjacency_.begin(), facets_adjacency_.end(), [](const auto& a, const auto& b) {
        if (a[0] != b[0]) {
            return a[0] < b[0];
        }
        return a[1] < b[1];
    });
    facets_adjacency_.erase(std::unique(facets_adjacency_.begin(), facets_adjacency_.end()), facets_adjacency_.end());
}

void Mesh::compute_interface_faces() {
    interfaces_.clear();
    interfacets_.clear();

    std::vector<int> interfacet_boundary_edges;
    for (size_t ei = 0; ei < edges_faces_.size(); ++ei) {
        if (edges_faces_[ei].size() > 2) {
            interfacet_boundary_edges.push_back(static_cast<int>(ei));
        }
    }

    for (size_t fi = 0; fi < facets_boundary_edges_.size(); ++fi) {
        const auto& b = facets_boundary_edges_[fi];
        if (b.empty()) {
            continue;
        }
        bool all_internal_boundary = true;
        for (const int e : b) {
            if (std::find(interfacet_boundary_edges.begin(), interfacet_boundary_edges.end(), e) == interfacet_boundary_edges.end()) {
                all_internal_boundary = false;
                break;
            }
        }
        if (all_internal_boundary) {
            interfacets_.push_back(static_cast<int>(fi));
            for (const int face : facets_[fi].faces) {
                interfaces_.push_back(face);
            }
        }
    }
    std::sort(interfaces_.begin(), interfaces_.end());
    interfaces_.erase(std::unique(interfaces_.begin(), interfaces_.end()), interfaces_.end());
}

void Mesh::compute_enclosed_volume() {
    const Vec3 center = mul(add(bounds_min_, bounds_max_), 0.5);
    double v = 0.0;
    for (const auto& f : faces_) {
        const Vec3 a = sub(vertices_[f[0]], center);
        const Vec3 b = sub(vertices_[f[1]], center);
        const Vec3 c = sub(vertices_[f[2]], center);
        v += std::abs(dot(a, cross(b, c))) / 6.0;
    }
    volume_ = v;
}

void Mesh::ensure_volume_tetrahedra() const {
    if (simplices_ready_) {
        return;
    }
    const_cast<Mesh*>(this)->build_volume_tetrahedra();
}

void Mesh::build_volume_tetrahedra() {
    simplices_.clear();
    auto build_star_fallback = [this]() {
        Vec3 center {0.0, 0.0, 0.0};
        for (const auto& v : vertices_) {
            center = add(center, v);
        }
        center = mul(center, 1.0 / static_cast<double>(vertices_.size()));
        if (!contains_point_naive(center)) {
            center = mul(add(bounds_min_, bounds_max_), 0.5);
        }
        for (const auto& f : faces_) {
            const Vec3 a = vertices_[f[0]];
            const Vec3 b = vertices_[f[1]];
            const Vec3 c = vertices_[f[2]];
            const double vol = std::abs(dot(sub(a, center), cross(sub(b, center), sub(c, center)))) / 6.0;
            if (vol <= 1e-16) {
                continue;
            }
            Simplex s;
            s.v = {center, a, b, c};
            s.volume = vol;
            s.bounds = {min_vec(min_vec(center, a), min_vec(b, c)), max_vec(max_vec(center, a), max_vec(b, c))};
            const std::array<std::array<double, 3>, 3> A {{
                {{a[0] - center[0], b[0] - center[0], c[0] - center[0]}},
                {{a[1] - center[1], b[1] - center[1], c[1] - center[1]}},
                {{a[2] - center[2], b[2] - center[2], c[2] - center[2]}}
            }};
            bool ok = false;
            s.invA = inverse3x3(A, ok);
            if (!ok) {
                continue;
            }
            simplices_.push_back(s);
        }
    };

    const Vec3 ext = sub(bounds_max_, bounds_min_);
    const double min_extent = std::max(1e-9, std::min({ext[0], ext[1], ext[2]}));
    const double max_edge_div = std::pow(10.0, std::ceil(std::log10(min_extent)) - 1.0);
    const double dedup_h = std::max(1e-10, max_edge_div * 1e-3);

    std::vector<Vec3> points = vertices_;
    points.reserve(vertices_.size() + edges_.size() * 4);

    for (const auto& e : edges_) {
        const Vec3 a = vertices_[e[0]];
        const Vec3 b = vertices_[e[1]];
        const Vec3 ab = sub(b, a);
        const double len = norm(ab);
        const int n = std::max(1, static_cast<int>(std::ceil(len / max_edge_div)));
        for (int j = 1; j < n; ++j) {
            points.push_back(add(a, mul(ab, static_cast<double>(j) / static_cast<double>(n))));
        }
    }

    const int nx = std::max(1, static_cast<int>(std::ceil(ext[0] / max_edge_div)));
    const int ny = std::max(1, static_cast<int>(std::ceil(ext[1] / max_edge_div)));
    const int nz = std::max(1, static_cast<int>(std::ceil(ext[2] / max_edge_div)));
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                const Vec3 p {
                    bounds_min_[0] + (static_cast<double>(ix) + 0.5) * ext[0] / static_cast<double>(nx),
                    bounds_min_[1] + (static_cast<double>(iy) + 0.5) * ext[1] / static_cast<double>(ny),
                    bounds_min_[2] + (static_cast<double>(iz) + 0.5) * ext[2] / static_cast<double>(nz)
                };
                if (!contains_point_naive(p)) {
                    continue;
                }
                const auto cp = nearest_point(p);
                if (cp.distance > max_edge_div * 0.5) {
                    points.push_back(p);
                }
            }
        }
    }

    // Surface samples by ray casting from outside, similar to Python sample_surface_points path.
    for (int d = 0; d < 3; ++d) {
        const int a = (d + 1) % 3;
        const int b = (d + 2) % 3;
        const int na = std::max(1, static_cast<int>(std::ceil(ext[a] / max_edge_div)));
        const int nb = std::max(1, static_cast<int>(std::ceil(ext[b] / max_edge_div)));
        std::vector<Vec3> ray_o;
        std::vector<Vec3> ray_v;
        ray_o.reserve(static_cast<size_t>(na * nb));
        ray_v.reserve(static_cast<size_t>(na * nb));
        for (int ia = 0; ia < na; ++ia) {
            for (int ib = 0; ib < nb; ++ib) {
                Vec3 p {0.0, 0.0, 0.0};
                p[d] = bounds_min_[d] - 1e-8 * std::max(1.0, ext[d]);
                p[a] = bounds_min_[a] + (static_cast<double>(ia) + 0.5) * ext[a] / static_cast<double>(na);
                p[b] = bounds_min_[b] + (static_cast<double>(ib) + 0.5) * ext[b] / static_cast<double>(nb);
                Vec3 v {0.0, 0.0, 0.0};
                v[d] = 1.0;
                ray_o.push_back(p);
                ray_v.push_back(v);
            }
        }
        const auto [xc, _, __, ___] = trace_boundary_intersections(ray_o, ray_v);
        points.insert(points.end(), xc.begin(), xc.end());
    }

    std::unordered_map<std::string, int> uniq_map;
    std::vector<Vec3> uniq_points;
    uniq_points.reserve(points.size());
    for (const auto& p : points) {
        const auto key = quant_key(p, dedup_h);
        if (uniq_map.emplace(key, static_cast<int>(uniq_points.size())).second) {
            uniq_points.push_back(p);
        }
    }
    if (uniq_points.size() < 4) {
        build_star_fallback();
        simplices_ready_ = true;
        return;
    }

    namespace fs = std::filesystem;
    std::mt19937_64 tmp_rng(std::random_device{}());
    const auto tmp_tag = std::to_string(static_cast<unsigned long long>(tmp_rng()));
    const fs::path tmp = fs::temp_directory_path() / ("epmc_qhull_" + tmp_tag + ".txt");
    {
        std::ofstream in(tmp);
        in << "3 " << uniq_points.size() << '\n';
        for (const auto& p : uniq_points) {
            in << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
        }
    }

    std::string output;
    const std::string cmd = "qdelaunay QJ i < \"" + tmp.string() + "\"";
    FILE* pipe = popen(cmd.c_str(), "r");
    int status = -1;
    if (pipe != nullptr) {
        char buffer[4096];
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            output += buffer;
        }
        status = pclose(pipe);
    }
    std::error_code ec;
    fs::remove(tmp, ec);

    if (status != 0 || output.empty()) {
        build_star_fallback();
        simplices_ready_ = true;
        return;
    }

    std::istringstream iss(output);
    int ns = 0;
    iss >> ns;
    if (!iss || ns <= 0) {
        build_star_fallback();
        simplices_ready_ = true;
        return;
    }

    for (int k = 0; k < ns; ++k) {
        int i0 = -1, i1 = -1, i2 = -1, i3 = -1;
        iss >> i0 >> i1 >> i2 >> i3;
        if (!iss) {
            break;
        }
        if (i0 < 0 || i1 < 0 || i2 < 0 || i3 < 0 ||
            i0 >= static_cast<int>(uniq_points.size()) ||
            i1 >= static_cast<int>(uniq_points.size()) ||
            i2 >= static_cast<int>(uniq_points.size()) ||
            i3 >= static_cast<int>(uniq_points.size())) {
            continue;
        }
        const Vec3 a = uniq_points[i0];
        const Vec3 b = uniq_points[i1];
        const Vec3 c = uniq_points[i2];
        const Vec3 d = uniq_points[i3];
        const Vec3 centroid = mul(add(add(a, b), add(c, d)), 0.25);
        if (!contains_point_naive(centroid)) {
            continue;
        }
        const double vol = std::abs(dot(sub(a, d), cross(sub(b, d), sub(c, d)))) / 6.0;
        if (vol <= 1e-9) {
            continue;
        }

        Simplex s;
        s.v = {a, b, c, d};
        s.volume = vol;
        s.bounds = {min_vec(min_vec(a, b), min_vec(c, d)), max_vec(max_vec(a, b), max_vec(c, d))};
        const std::array<std::array<double, 3>, 3> A {{
            {{b[0] - a[0], c[0] - a[0], d[0] - a[0]}},
            {{b[1] - a[1], c[1] - a[1], d[1] - a[1]}},
            {{b[2] - a[2], c[2] - a[2], d[2] - a[2]}}
        }};
        bool ok = false;
        s.invA = inverse3x3(A, ok);
        if (!ok) {
            continue;
        }
        simplices_.push_back(s);
    }

    if (simplices_.empty()) {
        build_star_fallback();
    }
    simplices_ready_ = true;
}

int Mesh::nearest_facet(const Vec3& p) const {
    const auto fq = nearest_face(p);
    if (fq.index < 0) {
        return -1;
    }
    return face_to_facet_[fq.index];
}

Mesh::FaceQuery Mesh::nearest_face(const Vec3& p) const {
    FaceQuery q;
    q.index = -1;
    q.distance = std::numeric_limits<double>::infinity();
    q.closest = p;
    for (size_t i = 0; i < faces_.size(); ++i) {
        const auto& f = faces_[i];
        const Vec3 cp = closest_point_on_triangle(p, vertices_[f[0]], vertices_[f[1]], vertices_[f[2]]);
        const double d = norm(sub(p, cp));
        if (d < q.distance) {
            q.distance = d;
            q.closest = cp;
            q.index = static_cast<int>(i);
        }
    }
    return q;
}

std::vector<Mesh::FaceQuery> Mesh::nearest_face(const std::vector<Vec3>& p) const {
    std::vector<FaceQuery> out;
    out.reserve(p.size());
    for (const auto& pt : p) {
        out.push_back(nearest_face(pt));
    }
    return out;
}

Mesh::EdgeQuery Mesh::nearest_edge(const Vec3& p) const {
    EdgeQuery q;
    q.index = -1;
    q.distance = std::numeric_limits<double>::infinity();
    q.closest = p;
    for (size_t i = 0; i < edges_.size(); ++i) {
        const auto& e = edges_[i];
        const Vec3 a = vertices_[e[0]];
        const Vec3 b = vertices_[e[1]];
        const Vec3 ab = sub(b, a);
        const double denom = dot(ab, ab);
        double t = 0.0;
        if (denom > 0.0) {
            t = dot(sub(p, a), ab) / denom;
        }
        t = std::clamp(t, 0.0, 1.0);
        const Vec3 cp = add(a, mul(ab, t));
        const double d = norm(sub(p, cp));
        if (d < q.distance) {
            q.distance = d;
            q.closest = cp;
            q.index = static_cast<int>(i);
        }
    }
    return q;
}

std::vector<Mesh::EdgeQuery> Mesh::nearest_edge(const std::vector<Vec3>& p) const {
    std::vector<EdgeQuery> out;
    out.reserve(p.size());
    for (const auto& pt : p) {
        out.push_back(nearest_edge(pt));
    }
    return out;
}

Mesh::PointQuery Mesh::nearest_point(const Vec3& p) const {
    const auto fq = nearest_face(p);
    const auto eq = nearest_edge(p);
    PointQuery out;
    if (eq.distance < fq.distance) {
        out.closest = eq.closest;
        out.distance = eq.distance;
        out.from_edge = true;
    } else {
        out.closest = fq.closest;
        out.distance = fq.distance;
        out.from_edge = false;
    }
    return out;
}

std::vector<Mesh::PointQuery> Mesh::nearest_point(const std::vector<Vec3>& p) const {
    std::vector<PointQuery> out;
    out.reserve(p.size());
    for (const auto& pt : p) {
        out.push_back(nearest_point(pt));
    }
    return out;
}

bool Mesh::contains_point_naive(const Vec3& p) const {
    // Ray casting in +x direction.
    const Vec3 dir {1.0, 0.0, 0.0};
    int hits = 0;
    for (const auto& f : faces_) {
        double t = 0.0;
        if (ray_intersects_triangle(p, dir, vertices_[f[0]], vertices_[f[1]], vertices_[f[2]], t)) {
            ++hits;
        }
    }
    return (hits % 2) == 1;
}

bool Mesh::contains_point(const Vec3& p) const {
    ensure_volume_tetrahedra();
    if (simplices_.empty()) {
        return contains_point_naive(p);
    }
    constexpr double tol = 1e-10;
    for (const auto& s : simplices_) {
        if (p[0] < s.bounds[0][0] - tol || p[1] < s.bounds[0][1] - tol || p[2] < s.bounds[0][2] - tol ||
            p[0] > s.bounds[1][0] + tol || p[1] > s.bounds[1][1] + tol || p[2] > s.bounds[1][2] + tol) {
            continue;
        }
        const Vec3 b = matvec3(s.invA, sub(p, s.v[0]));
        const double l1 = b[0];
        const double l2 = b[1];
        const double l3 = b[2];
        const double l0 = 1.0 - l1 - l2 - l3;
        if (l0 >= -tol && l1 >= -tol && l2 >= -tol && l3 >= -tol &&
            l0 <= 1.0 + tol && l1 <= 1.0 + tol && l2 <= 1.0 + tol && l3 <= 1.0 + tol) {
            return true;
        }
    }
    return false;
}

std::vector<bool> Mesh::contains_point(const std::vector<Vec3>& p) const {
    std::vector<bool> out;
    out.reserve(p.size());
    for (const auto& pt : p) {
        out.push_back(contains_point(pt));
    }
    return out;
}

std::tuple<Vec3, double, int> Mesh::trace_boundary_intersection(const Vec3& x, const Vec3& v) const {
    constexpr double eps = 1e-10;
    double best_t = std::numeric_limits<double>::infinity();
    int best_face = -1;
    Vec3 best_p {0.0, 0.0, 0.0};

    for (size_t i = 0; i < faces_.size(); ++i) {
        const auto denom = dot(v, face_normals_[i]);
        if (std::abs(denom) < eps) {
            continue;
        }
        const double t = -(dot(x, face_normals_[i]) + face_k_[i]) / denom;
        if (t <= eps || t >= best_t) {
            continue;
        }
        const Vec3 p = add(x, mul(v, t));
        const auto& bb = face_bounds_[i];
        if (p[0] < bb[0][0] - 1e-9 || p[1] < bb[0][1] - 1e-9 || p[2] < bb[0][2] - 1e-9 ||
            p[0] > bb[1][0] + 1e-9 || p[1] > bb[1][1] + 1e-9 || p[2] > bb[1][2] + 1e-9) {
            continue;
        }
        const auto& f = faces_[i];
        const Vec3 cp = closest_point_on_triangle(p, vertices_[f[0]], vertices_[f[1]], vertices_[f[2]]);
        if (norm(sub(cp, p)) > 1e-7) {
            continue;
        }
        best_t = t;
        best_p = p;
        best_face = static_cast<int>(i);
    }

    if (best_face < 0) {
        return {x, std::numeric_limits<double>::infinity(), -1};
    }
    return {best_p, best_t, face_to_facet_[best_face]};
}

std::tuple<std::vector<Vec3>, std::vector<double>, std::vector<int>, std::vector<int>> Mesh::trace_boundary_intersections(
    const std::vector<Vec3>& x,
    const std::vector<Vec3>& v) const {
    if (x.size() != v.size()) {
        throw std::runtime_error("trace_boundary_intersections requires x.size()==v.size().");
    }
    constexpr double tol = 1e-10;
    std::vector<Vec3> current = x;
    std::vector<bool> active(x.size(), true);
    std::vector<Vec3> xc;
    std::vector<double> tc;
    std::vector<int> fc;
    std::vector<int> ic;
    xc.reserve(x.size());
    tc.reserve(x.size());
    fc.reserve(x.size());
    ic.reserve(x.size());

    bool any_active = true;
    while (any_active) {
        any_active = false;
        for (size_t i = 0; i < x.size(); ++i) {
            if (!active[i]) {
                continue;
            }
            any_active = true;
            const auto [p, t, f] = trace_boundary_intersection(current[i], v[i]);
            if (f >= 0 && !std::isinf(t)) {
                xc.push_back(p);
                tc.push_back(t);
                fc.push_back(f);
                ic.push_back(static_cast<int>(i));
                current[i] = add(p, mul(v[i], tol));
            } else {
                active[i] = false;
            }
        }
    }
    return {xc, tc, fc, ic};
}

std::vector<Vec3> Mesh::sample_surface_points(int n, const std::vector<int>& facets) const {
    if (n <= 0) {
        return {};
    }
    std::vector<int> face_pool;
    if (facets.empty()) {
        face_pool.resize(faces_.size());
        for (size_t i = 0; i < faces_.size(); ++i) {
            face_pool[i] = static_cast<int>(i);
        }
    } else {
        for (const int fct : facets) {
            if (fct < 0 || fct >= static_cast<int>(facets_.size())) {
                continue;
            }
            for (const int fi : facets_[fct].faces) {
                face_pool.push_back(fi);
            }
        }
    }
    if (face_pool.empty()) {
        throw std::runtime_error("No faces available for surface sampling.");
    }

    std::vector<double> w;
    w.reserve(face_pool.size());
    for (const int fi : face_pool) {
        w.push_back(std::max(face_areas_[fi], 0.0));
    }
    std::discrete_distribution<int> face_pick(w.begin(), w.end());
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    std::vector<Vec3> out;
    out.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const int pick = face_pool[face_pick(rng)];
        const auto& tri = faces_[pick];
        const Vec3 a = vertices_[tri[0]];
        const Vec3 b = vertices_[tri[1]];
        const Vec3 c = vertices_[tri[2]];
        const double r1 = uni(rng);
        const double r2 = uni(rng);
        const double sr1 = std::sqrt(r1);
        const double w0 = 1.0 - sr1;
        const double w1 = sr1 * (1.0 - r2);
        const double w2 = sr1 * r2;
        out.push_back({
            w0 * a[0] + w1 * b[0] + w2 * c[0],
            w0 * a[1] + w1 * b[1] + w2 * c[1],
            w0 * a[2] + w1 * b[2] + w2 * c[2]
        });
    }
    return out;
}

std::vector<Vec3> Mesh::sample_volume_points_naive(int n) const {
    if (n <= 0) {
        return {};
    }
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> ux(bounds_min_[0], bounds_max_[0]);
    std::uniform_real_distribution<double> uy(bounds_min_[1], bounds_max_[1]);
    std::uniform_real_distribution<double> uz(bounds_min_[2], bounds_max_[2]);

    std::vector<Vec3> out;
    out.reserve(static_cast<size_t>(n));
    const int max_trials = std::max(10000, n * 200);
    int trials = 0;
    while (static_cast<int>(out.size()) < n && trials < max_trials) {
        ++trials;
        const Vec3 p {ux(rng), uy(rng), uz(rng)};
        if (contains_point_naive(p)) {
            out.push_back(p);
        }
    }
    if (static_cast<int>(out.size()) < n) {
        throw std::runtime_error("sample_volume_points_naive: insufficient accepted samples.");
    }
    return out;
}

std::vector<Vec3> Mesh::sample_volume_points(int n) const {
    if (n <= 0) {
        return {};
    }
    ensure_volume_tetrahedra();
    if (simplices_.empty()) {
        return sample_volume_points_naive(n);
    }

    std::vector<double> w;
    w.reserve(simplices_.size());
    for (const auto& s : simplices_) {
        w.push_back(std::max(s.volume, 0.0));
    }
    const double wsum = std::accumulate(w.begin(), w.end(), 0.0);
    if (wsum <= 0.0) {
        return sample_volume_points_naive(n);
    }

    std::discrete_distribution<int> simplex_pick(w.begin(), w.end());
    std::mt19937_64 rng(std::random_device{}());
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    std::vector<Vec3> out;
    out.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const auto& s = simplices_[simplex_pick(rng)];
        // Dirichlet(1,1,1,1) via normalized exponential.
        double a0 = -std::log(std::max(uni(rng), 1e-15));
        double a1 = -std::log(std::max(uni(rng), 1e-15));
        double a2 = -std::log(std::max(uni(rng), 1e-15));
        double a3 = -std::log(std::max(uni(rng), 1e-15));
        const double sum = a0 + a1 + a2 + a3;
        a0 /= sum;
        a1 /= sum;
        a2 /= sum;
        a3 /= sum;
        out.push_back({
            a0 * s.v[0][0] + a1 * s.v[1][0] + a2 * s.v[2][0] + a3 * s.v[3][0],
            a0 * s.v[0][1] + a1 * s.v[1][1] + a2 * s.v[2][1] + a3 * s.v[3][1],
            a0 * s.v[0][2] + a1 * s.v[1][2] + a2 * s.v[2][2] + a3 * s.v[3][2]
        });
    }
    return out;
}

int Mesh::facet_index_for_face(int face_index) const {
    if (face_index < 0 || face_index >= static_cast<int>(face_to_facet_.size())) {
        return -1;
    }
    return face_to_facet_[face_index];
}

std::vector<int> Mesh::facet_indices_for_faces(const std::vector<int>& face_indices) const {
    std::vector<int> out;
    out.reserve(face_indices.size());
    for (const int idx : face_indices) {
        out.push_back(facet_index_for_face(idx));
    }
    return out;
}

void Mesh::erase_vertices(const std::vector<int>& indices) {
    if (indices.empty()) {
        return;
    }
    std::set<int> rm(indices.begin(), indices.end());
    std::vector<int> old_to_new(vertices_.size(), -1);
    std::vector<Vec3> new_vertices;
    new_vertices.reserve(vertices_.size() - std::min(vertices_.size(), rm.size()));
    for (size_t i = 0; i < vertices_.size(); ++i) {
        if (rm.count(static_cast<int>(i)) == 0) {
            old_to_new[i] = static_cast<int>(new_vertices.size());
            new_vertices.push_back(vertices_[i]);
        }
    }
    std::vector<Tri> new_faces;
    new_faces.reserve(faces_.size());
    for (const auto& f : faces_) {
        if (old_to_new[f[0]] < 0 || old_to_new[f[1]] < 0 || old_to_new[f[2]] < 0) {
            continue;
        }
        new_faces.push_back({old_to_new[f[0]], old_to_new[f[1]], old_to_new[f[2]]});
    }
    if (new_vertices.empty() || new_faces.empty()) {
        throw std::runtime_error("erase_vertices produced empty mesh.");
    }
    vertices_ = std::move(new_vertices);
    faces_ = std::move(new_faces);
    rebuild_cached_properties();
}

void Mesh::erase_faces(const std::vector<int>& indices) {
    if (indices.empty()) {
        return;
    }
    std::set<int> rm(indices.begin(), indices.end());
    std::vector<Tri> new_faces;
    new_faces.reserve(faces_.size() - std::min(faces_.size(), rm.size()));
    for (size_t i = 0; i < faces_.size(); ++i) {
        if (rm.count(static_cast<int>(i)) == 0) {
            new_faces.push_back(faces_[i]);
        }
    }
    if (new_faces.empty()) {
        throw std::runtime_error("erase_faces produced empty mesh faces.");
    }
    faces_ = std::move(new_faces);
    rebuild_cached_properties();
}

void Mesh::remove_unreferenced_vertices() {
    std::vector<bool> used(vertices_.size(), false);
    for (const auto& f : faces_) {
        used[f[0]] = true;
        used[f[1]] = true;
        used[f[2]] = true;
    }
    std::vector<int> rm;
    for (size_t i = 0; i < used.size(); ++i) {
        if (!used[i]) {
            rm.push_back(static_cast<int>(i));
        }
    }
    if (!rm.empty()) {
        erase_vertices(rm);
    }
}

void Mesh::merge_duplicate_vertices(double tol) {
    if (vertices_.empty() || faces_.empty()) {
        return;
    }
    // Deduplicate vertices by quantization.
    auto q = [tol](double x) -> std::int64_t {
        return static_cast<std::int64_t>(std::llround(x / tol));
    };
    std::unordered_map<std::string, int> vmap;
    std::vector<int> old_to_new(vertices_.size(), -1);
    std::vector<Vec3> new_vertices;
    new_vertices.reserve(vertices_.size());
    for (size_t i = 0; i < vertices_.size(); ++i) {
        const auto& v = vertices_[i];
        const std::string key = std::to_string(q(v[0])) + ":" + std::to_string(q(v[1])) + ":" + std::to_string(q(v[2]));
        auto it = vmap.find(key);
        if (it == vmap.end()) {
            const int idx = static_cast<int>(new_vertices.size());
            vmap.emplace(key, idx);
            old_to_new[i] = idx;
            new_vertices.push_back(v);
        } else {
            old_to_new[i] = it->second;
        }
    }
    std::vector<Tri> reindexed_faces;
    reindexed_faces.reserve(faces_.size());
    for (const auto& f : faces_) {
        Tri r {old_to_new[f[0]], old_to_new[f[1]], old_to_new[f[2]]};
        if (r[0] == r[1] || r[0] == r[2] || r[1] == r[2]) {
            continue;
        }
        reindexed_faces.push_back(r);
    }
    // Deduplicate faces ignoring winding.
    std::set<std::array<int, 3>> uniq;
    std::vector<Tri> dedup_faces;
    dedup_faces.reserve(reindexed_faces.size());
    for (const auto& f : reindexed_faces) {
        std::array<int, 3> s = f;
        std::sort(s.begin(), s.end());
        if (uniq.insert(s).second) {
            dedup_faces.push_back(f);
        }
    }
    if (new_vertices.empty() || dedup_faces.empty()) {
        throw std::runtime_error("merge_duplicate_vertices produced empty mesh.");
    }
    vertices_ = std::move(new_vertices);
    faces_ = std::move(dedup_faces);
    remove_unreferenced_vertices();
    rebuild_cached_properties();
}

void Mesh::write_stl(const std::string& path, const std::string& name) const {
    namespace fs = std::filesystem;
    fs::create_directories(path);
    const fs::path out_path = fs::path(path) / (name + ".stl");
    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("Failed to open STL output: " + out_path.string());
    }
    out << "solid " << name << '\n';
    out.setf(std::ios::scientific);
    out.precision(6);
    for (size_t i = 0; i < faces_.size(); ++i) {
        const auto& f = faces_[i];
        out << "facet normal " << face_normals_[i][0] << " " << face_normals_[i][1] << " " << face_normals_[i][2] << '\n';
        out << "    outer loop\n";
        out << "        vertex " << vertices_[f[0]][0] << " " << vertices_[f[0]][1] << " " << vertices_[f[0]][2] << '\n';
        out << "        vertex " << vertices_[f[1]][0] << " " << vertices_[f[1]][1] << " " << vertices_[f[1]][2] << '\n';
        out << "        vertex " << vertices_[f[2]][0] << " " << vertices_[f[2]][1] << " " << vertices_[f[2]][2] << '\n';
        out << "    endloop\n";
        out << "endfacet\n";
    }
    out << "endsolid " << name << '\n';
}
