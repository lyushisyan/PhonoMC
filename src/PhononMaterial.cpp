#include "PhononMaterial.h"

#include "material/BzQuadrature.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace {
using Mat3 = std::array<std::array<double, 3>, 3>;
using Vec3 = std::array<double, 3>;

// 函数说明：将分数 q 点映射到倒易空间 k 向量。
Vec3 q_to_k(const Mat3& reciprocal_lattice, const Vec3& q) {
    return {reciprocal_lattice[0][0] * q[0] + reciprocal_lattice[0][1] * q[1] +
                reciprocal_lattice[0][2] * q[2],
            reciprocal_lattice[1][0] * q[0] + reciprocal_lattice[1][1] * q[1] +
                reciprocal_lattice[1][2] * q[2],
            reciprocal_lattice[2][0] * q[0] + reciprocal_lattice[2][1] * q[1] +
                reciprocal_lattice[2][2] * q[2]};
}

// 函数说明：计算三维向量模长。
double vec_norm(const Vec3& v) {
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}
} // namespace

// File-free material model. Application loading policy lives in MaterialFactory.
PhononMaterial::PhononMaterial(phonomc::MaterialData data,
                               double temperature_lookup_dt,
                               bool include_stationary_modes,
                               bool complete_callaway_bz)
    : temperature_lookup_dt_(temperature_lookup_dt),
      include_stationary_modes_(include_stationary_modes) {
    if (!std::isfinite(temperature_lookup_dt_) || temperature_lookup_dt_ <= 0.0) {
        throw std::invalid_argument("Temperature lookup spacing must be finite and positive");
    }
    phonomc::validate_material_data(data);
    if (complete_callaway_bz && !data.gamma_normal.empty()) {
        bz_refinement_ = phonomc::complete_callaway_bz(data);
    }
    initialize_material_data(std::move(data));
    initialize_temperature_lookup();
}

PhononMaterial::PhononMaterial(SyntheticTag, double temperature_lookup_dt)
    : temperature_lookup_dt_(temperature_lookup_dt) {
    if (!std::isfinite(temperature_lookup_dt_) || temperature_lookup_dt_ <= 0.0) {
        throw std::invalid_argument("Temperature lookup spacing must be finite and positive");
    }
    build_synthetic_test_modes();
    initialize_temperature_lookup();
}

PhononMaterial PhononMaterial::synthetic_for_testing(double temperature_lookup_dt) {
    return PhononMaterial(SyntheticTag{}, temperature_lookup_dt);
}

int PhononMaterial::flatten_mode_index(const Mode& mode) const {
    return mode[0] * branch_count_ + mode[1];
}

// 函数说明：按活跃模态索引返回具体 (q, branch) 模式。
PhononMaterial::Mode PhononMaterial::active_mode_at(int active_index) const {
    if (active_index < 0 || active_index >= static_cast<int>(active_mode_list_.size())) {
        return {0, 0};
    }
    return active_mode_list_[static_cast<size_t>(active_index)];
}

// 函数说明：将模态反查为活跃模态索引。
int PhononMaterial::active_index_for_mode(const Mode& mode) const {
    const int fi = flatten_mode_index(mode);
    if (fi < 0 || fi >= static_cast<int>(flat_to_active_index_.size())) {
        return -1;
    }
    return flat_to_active_index_[static_cast<size_t>(fi)];
}

// 函数说明：查询简并模态对应的伴随支。
int PhononMaterial::degenerate_partner_branch(const Mode& mode) const {
    const int fi = flatten_mode_index(mode);
    if (fi < 0 || fi >= static_cast<int>(degenerate_partner_branch_data_.size())) {
        return -1;
    }
    return degenerate_partner_branch_data_[static_cast<size_t>(fi)];
}

void PhononMaterial::initialize_material_data(phonomc::MaterialData data) {
    unit_cell_volume_ = data.lattice.volume_a3;
    qpoint_count_ = data.qpoints;
    branch_count_ = data.branches;
    const int nm = qpoint_count_ * branch_count_;
    mode_angular_frequency_data_.resize(static_cast<size_t>(nm));
    mode_wavevector_norm_data_.assign(static_cast<size_t>(nm), 0.0);
    mode_angular_frequency_window_data_.assign(static_cast<size_t>(nm), 0.0);
    for (int i = 0; i < nm; ++i) {
        mode_angular_frequency_data_[static_cast<size_t>(i)] =
            std::max(0.0, data.frequency_thz[static_cast<size_t>(i)]) * 2.0 * M_PI;
    }
    mode_group_velocity_data_ = std::move(data.group_velocity);
    temperature_samples_ = std::move(data.temperatures);
    gamma_table_ = std::move(data.gamma);
    gamma_normal_ = std::move(data.gamma_normal);
    gamma_umklapp_ = std::move(data.gamma_umklapp);
    gamma_extra_.assign(gamma_table_.size(), 0.0);
    for (const auto& [name, values] : data.gamma_resistive) {
        for (size_t i = 0; i < gamma_extra_.size(); ++i) {
            gamma_extra_[i] += values[values.size() == static_cast<size_t>(nm) ? i % nm : i];
        }
    }
    gamma_total_ = gamma_table_;
    for (size_t i = 0; i < gamma_total_.size(); ++i) {
        gamma_total_[i] += gamma_extra_[i];
    }
    qpoint_fractions_ = data.qpoint_fractions;
    reciprocal_ = data.lattice.reciprocal;
    mesh_ = data.mesh;
    mode_wavevector_data_.resize(static_cast<size_t>(nm));
    const Vec3 mesh_q{1.0 / (2.0 * static_cast<double>(data.mesh[0])),
                      1.0 / (2.0 * static_cast<double>(data.mesh[1])),
                      1.0 / (2.0 * static_cast<double>(data.mesh[2]))};
    const Vec3 k_grid = q_to_k(data.lattice.reciprocal, mesh_q);
    for (int q = 0; q < qpoint_count_; ++q) {
        const Vec3 kvec =
            q_to_k(data.lattice.reciprocal, data.qpoint_fractions[static_cast<size_t>(q)]);
        const double kn = vec_norm(kvec);
        for (int b = 0; b < branch_count_; ++b) {
            const int fi = q * branch_count_ + b;
            const Vec3 v = mode_group_velocity_data_[static_cast<size_t>(fi)];
            mode_wavevector_norm_data_[static_cast<size_t>(fi)] = kn;
            mode_wavevector_data_[static_cast<size_t>(fi)] = kvec;
            mode_angular_frequency_window_data_[static_cast<size_t>(fi)] = std::sqrt(
                (v[0] * k_grid[0]) * (v[0] * k_grid[0]) + (v[1] * k_grid[1]) * (v[1] * k_grid[1]) +
                (v[2] * k_grid[2]) * (v[2] * k_grid[2]));
        }
    }

    active_mode_list_.clear();
    // Keep replicas adjacent in the sampled active bank. The q-major raw tables
    // remain unchanged. This also aligns the underlying physical-mode ordering
    // with the unrefined bank for paired random-stream comparisons.
    for (int base_q = 0; base_q < qpoint_count_; base_q += bz_refinement_) {
        for (int b = 0; b < branch_count_; ++b) {
            for (unsigned copy = 0; copy < bz_refinement_; ++copy) {
                const int q = base_q + copy;
                const auto& v =
                    mode_group_velocity_data_[static_cast<size_t>(q * branch_count_ + b)];
                const double vn = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                if (include_stationary_modes_
                        ? mode_angular_frequency_data_[static_cast<size_t>(q * branch_count_ + b)] >
                              0.0
                        : vn > 1e-12) {
                    active_mode_list_.push_back({q, b});
                }
            }
        }
    }
    active_mode_count_ = static_cast<int>(active_mode_list_.size());
    if (active_mode_count_ == 0) {
        throw std::runtime_error("No active modes found in material data.");
    }
    build_active_mode_maps();
}

// 函数说明：构建测试用合成模态库（非物理生产用途）。
void PhononMaterial::build_synthetic_test_modes() {
    const int n_q = 64;
    qpoint_count_ = n_q;
    branch_count_ = 3;
    active_mode_count_ = n_q * branch_count_;
    unit_cell_volume_ = 1.0;
    mode_angular_frequency_data_.resize(static_cast<size_t>(active_mode_count_));
    mode_wavevector_norm_data_.assign(static_cast<size_t>(active_mode_count_), 1.0);
    mode_angular_frequency_window_data_.assign(static_cast<size_t>(active_mode_count_), 1e-3);
    mode_group_velocity_data_.resize(static_cast<size_t>(active_mode_count_));
    temperature_samples_ = {100.0, 200.0, 300.0, 400.0, 500.0};
    gamma_table_.assign(temperature_samples_.size() * static_cast<size_t>(active_mode_count_),
                        0.05);
    gamma_total_ = gamma_table_;
    active_mode_list_.clear();
    active_mode_list_.reserve(static_cast<size_t>(active_mode_count_));

    std::mt19937_64 rng(20250209);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    for (int q = 0; q < n_q; ++q) {
        for (int b = 0; b < branch_count_; ++b) {
            const int m = q * branch_count_ + b;
            const double frac = (static_cast<double>(q) + 0.5) / static_cast<double>(n_q);
            mode_angular_frequency_data_[static_cast<size_t>(m)] =
                (2.0 + 18.0 * frac) * (1.0 + 0.15 * b) * 2.0 * M_PI;
            const double vmag = 20.0 + 90.0 * U(rng);
            const auto dir = random_unit_vector(rng);
            mode_group_velocity_data_[static_cast<size_t>(m)] = {
                dir[0] * vmag, dir[1] * vmag, dir[2] * vmag};
            active_mode_list_.push_back({q, b});
        }
    }
    build_active_mode_maps();
}

void PhononMaterial::build_active_mode_maps() {
    const int nm = qpoint_count_ * branch_count_;
    flat_to_active_index_.assign(static_cast<size_t>(nm), -1);
    for (size_t ai = 0; ai < active_mode_list_.size(); ++ai) {
        const int fi = flatten_mode_index(active_mode_list_[ai]);
        flat_to_active_index_[static_cast<size_t>(fi)] = static_cast<int>(ai);
    }
    build_degenerate_mode_map();
}

// 函数说明：采样随机单位方向用于合成模式或方向扰动。
PhononMaterial::Vec3 PhononMaterial::random_unit_vector(std::mt19937_64& rng) {
    std::normal_distribution<double> N(0.0, 1.0);
    Vec3 v{N(rng), N(rng), N(rng)};
    const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (n <= 1e-14) {
        return {1.0, 0.0, 0.0};
    }
    return {v[0] / n, v[1] / n, v[2] / n};
}

// 函数说明：在活跃模态集合中采样一个声子模式。
PhononMaterial::Mode PhononMaterial::sample_active_mode(std::mt19937_64& rng) const {
    if (active_mode_list_.empty()) {
        throw std::runtime_error("No active modes in PhononMaterial.");
    }
    std::uniform_int_distribution<int> pick(0, static_cast<int>(active_mode_list_.size()) - 1);
    return active_mode_list_[static_cast<size_t>(pick(rng))];
}

// 函数说明：查询给定模态的群速度矢量。
PhononMaterial::Vec3 PhononMaterial::mode_group_velocity(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_group_velocity_data_.size())) {
        return {0.0, 0.0, 0.0};
    }
    return mode_group_velocity_data_[static_cast<size_t>(m)];
}

// 函数说明：查询给定模态的角频率。
double PhononMaterial::mode_angular_frequency(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_angular_frequency_data_.size())) {
        return 1.0;
    }
    return mode_angular_frequency_data_[static_cast<size_t>(m)];
}

// 函数说明：查询给定模态波矢模长。
double PhononMaterial::mode_wavevector_norm(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_wavevector_norm_data_.size())) {
        return 0.0;
    }
    return mode_wavevector_norm_data_[static_cast<size_t>(m)];
}

// 函数说明：查询给定模态对应频率窗宽。
double PhononMaterial::mode_frequency_window(const Mode& mode) const {
    const int m = flatten_mode_index(mode);
    if (m < 0 || m >= static_cast<int>(mode_angular_frequency_window_data_.size())) {
        return 0.0;
    }
    return mode_angular_frequency_window_data_[static_cast<size_t>(m)];
}

// 函数说明：预计算简并模态映射以支持粗糙边界散射匹配。
void PhononMaterial::build_degenerate_mode_map() {
    const int nm = qpoint_count_ * branch_count_;
    degenerate_partner_branch_data_.assign(static_cast<size_t>(std::max(0, nm)), -1);
    if (qpoint_count_ <= 0 || branch_count_ <= 1) {
        return;
    }
    constexpr double tol = 1e-10;
    for (int q = 0; q < qpoint_count_; ++q) {
        for (int j1 = 0; j1 < branch_count_; ++j1) {
            const int f1 = q * branch_count_ + j1;
            for (int j2 = j1 + 1; j2 < branch_count_; ++j2) {
                const int f2 = q * branch_count_ + j2;
                if (std::abs(mode_angular_frequency_data_[static_cast<size_t>(f1)] -
                             mode_angular_frequency_data_[static_cast<size_t>(f2)]) < tol) {
                    degenerate_partner_branch_data_[static_cast<size_t>(f1)] = j2;
                    degenerate_partner_branch_data_[static_cast<size_t>(f2)] = j1;
                    break;
                }
            }
        }
    }
}
