#include "geometry/MeshReader.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace phonomc {
namespace {
using Vec3 = std::array<double, 3>;
using Tri = std::array<int, 3>;

// 函数说明：将关键字统一为小写，避免输入大小写差异影响流程。
std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

// 函数说明：清理输入文本两端空白，保证配置与数据解析的稳健性。
std::string trim(const std::string& s) {
    const auto b = std::find_if_not(s.begin(), s.end(), [](unsigned char c) {
        return std::isspace(c);
    });
    if (b == s.end()) {
        return "";
    }
    const auto e = std::find_if_not(s.rbegin(), s.rend(), [](unsigned char c) {
                       return std::isspace(c);
                   }).base();
    return std::string(b, e);
}

// 函数说明：将几何坐标量化到容差网格，用于顶点去重键生成。
std::array<long long, 3> quant3(const Vec3& p, double h = 1e-10) {
    return {static_cast<long long>(std::llround(p[0] / h)),
            static_cast<long long>(std::llround(p[1] / h)),
            static_cast<long long>(std::llround(p[2] / h))};
}

// 函数说明：对顶点做量化去重并返回去重后的顶点索引。
int add_vertex_dedup(const Vec3& p,
                     std::vector<Vec3>& vertices,
                     std::unordered_map<std::string, int>& vmap,
                     double h = 1e-10) {
    const auto q = quant3(p, h);
    const std::string k =
        std::to_string(q[0]) + ":" + std::to_string(q[1]) + ":" + std::to_string(q[2]);
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
SurfaceMeshData load_obj_mesh(const std::filesystem::path& p) {
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
SurfaceMeshData load_ascii_stl_mesh(const std::filesystem::path& p) {
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
    T v{};
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return v;
}

// 函数说明：读取 Binary STL 并转换为内部网格表示。
SurfaceMeshData load_binary_stl_mesh(const std::filesystem::path& p) {
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
        (void)read_binary<float>(in);
        (void)read_binary<float>(in);
        (void)read_binary<float>(in);
        Vec3 p0{read_binary<float>(in), read_binary<float>(in), read_binary<float>(in)};
        Vec3 p1{read_binary<float>(in), read_binary<float>(in), read_binary<float>(in)};
        Vec3 p2{read_binary<float>(in), read_binary<float>(in), read_binary<float>(in)};
        (void)read_binary<uint16_t>(in);
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

} // namespace

// 函数说明：按扩展名调度网格加载器，统一支持 OBJ/STL。
SurfaceMeshData read_surface_mesh(const std::filesystem::path& p) {
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

} // namespace phonomc
