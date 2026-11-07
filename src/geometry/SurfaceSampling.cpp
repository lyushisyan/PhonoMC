#include "SurfaceMesh.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>

using Vec3 = SurfaceMesh::Vec3;

// Distribution objects stay local to each original call. Preserve pool ordering
// and the order of face/simplex selection followed by coordinate draws.
// 函数说明：按面面积加权在边界表面采样粒子注入点。
std::vector<Vec3> SurfaceMesh::sample_surface_points(int n, const std::vector<int>& facets) const {
    std::mt19937_64 rng(std::random_device{}());
    return sample_surface_points(n, facets, rng);
}

std::vector<Vec3> SurfaceMesh::sample_surface_points(int n,
                                                     const std::vector<int>& facets,
                                                     std::mt19937_64& rng) const {
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
        out.push_back({w0 * a[0] + w1 * b[0] + w2 * c[0],
                       w0 * a[1] + w1 * b[1] + w2 * c[1],
                       w0 * a[2] + w1 * b[2] + w2 * c[2]});
    }
    return out;
}

// 函数说明：用拒绝采样在体域内生成粒子初始位置。
std::vector<Vec3> SurfaceMesh::sample_volume_points_naive(int n) const {
    std::mt19937_64 rng(std::random_device{}());
    return sample_volume_points_naive(n, rng);
}

std::vector<Vec3> SurfaceMesh::sample_volume_points_naive(int n, std::mt19937_64& rng) const {
    if (n <= 0) {
        return {};
    }
    std::uniform_real_distribution<double> ux(bounds_min_[0], bounds_max_[0]);
    std::uniform_real_distribution<double> uy(bounds_min_[1], bounds_max_[1]);
    std::uniform_real_distribution<double> uz(bounds_min_[2], bounds_max_[2]);

    std::vector<Vec3> out;
    out.reserve(static_cast<size_t>(n));
    // Promote before multiplication: 40 million carriers need a trial budget
    // larger than INT_MAX, even though the accepted sample count fits in int.
    const std::int64_t max_trials = std::max<std::int64_t>(10000, std::int64_t{n} * 200);
    std::int64_t trials = 0;
    while (static_cast<int>(out.size()) < n && trials < max_trials) {
        ++trials;
        const Vec3 p{ux(rng), uy(rng), uz(rng)};
        if (contains_point_ray_cast(p)) {
            out.push_back(p);
        }
    }
    if (static_cast<int>(out.size()) < n) {
        throw std::runtime_error("sample_volume_points_naive: insufficient accepted samples.");
    }
    return out;
}

// 函数说明：优先基于四面体体积分布进行体内高效采样。
std::vector<Vec3> SurfaceMesh::sample_volume_points(int n) const {
    std::mt19937_64 rng(std::random_device{}());
    return sample_volume_points(n, rng);
}

std::vector<Vec3> SurfaceMesh::sample_volume_points(int n, std::mt19937_64& rng) const {
    if (n <= 0) {
        return {};
    }
    ensure_volume_tetrahedra();
    if (simplices_.empty()) {
        return sample_volume_points_naive(n, rng);
    }

    std::vector<double> w;
    w.reserve(simplices_.size());
    for (const auto& s : simplices_) {
        w.push_back(std::max(s.volume, 0.0));
    }
    const double wsum = std::accumulate(w.begin(), w.end(), 0.0);
    if (wsum <= 0.0) {
        return sample_volume_points_naive(n, rng);
    }

    std::discrete_distribution<int> simplex_pick(w.begin(), w.end());
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
        out.push_back({a0 * s.v[0][0] + a1 * s.v[1][0] + a2 * s.v[2][0] + a3 * s.v[3][0],
                       a0 * s.v[0][1] + a1 * s.v[1][1] + a2 * s.v[2][1] + a3 * s.v[3][1],
                       a0 * s.v[0][2] + a1 * s.v[1][2] + a2 * s.v[2][2] + a3 * s.v[3][2]});
    }
    return out;
}
