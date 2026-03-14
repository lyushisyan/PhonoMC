#include "MonteCarloSolver.h"

#include "SimulationDomain.h"
#include "PhononMaterial.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// 函数说明：执行三维向量加法，服务于粒子位置与速度相关更新。
MonteCarloSolver::Vec3 MonteCarloSolver::add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
// 函数说明：执行三维向量减法，服务于碰撞几何与距离计算。
MonteCarloSolver::Vec3 MonteCarloSolver::sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
// 函数说明：执行向量与标量缩放，用于时间推进与反射修正。
MonteCarloSolver::Vec3 MonteCarloSolver::mul(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
// 函数说明：计算向量点积，用于投影、入射角与法向分量判断。
double MonteCarloSolver::dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
// 函数说明：计算向量模长，用于速度归一化与阈值保护。
double MonteCarloSolver::norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// 函数说明：按线程数与网格数复用线程局部缓冲，降低 OpenMP 临时分配开销。
void MonteCarloSolver::ensure_tls_buffers(int thread_count, int nsv) {
    thread_count = std::max(1, thread_count);
    nsv = std::max(1, nsv);
    if (thread_count == tls_thread_count_ && nsv == tls_nsv_) {
        return;
    }
    tls_thread_count_ = thread_count;
    tls_nsv_ = nsv;
    const size_t total = static_cast<size_t>(thread_count) * static_cast<size_t>(nsv);
    energy_tls_buffer_.assign(total, 0.0);
    count_tls_buffer_.assign(total, 0);
    flux_tls_buffer_.assign(total, Vec3 {0.0, 0.0, 0.0});
}

// 函数说明：提供线程独立随机数发生器，保证并行采样过程的线程安全。
std::mt19937_64& MonteCarloSolver::thread_rng() const {
#ifdef _OPENMP
    static thread_local std::mt19937_64 tl_rng;
    static thread_local bool tl_seeded = false;
    if (!tl_seeded) {
        const unsigned tid = static_cast<unsigned>(omp_get_thread_num());
        std::seed_seq seq {
            static_cast<unsigned>(rng_seed_base_ & 0xffffffffu),
            static_cast<unsigned>((rng_seed_base_ >> 32) & 0xffffffffu),
            tid
        };
        tl_rng.seed(seq);
        tl_seeded = true;
    }
    return tl_rng;
#else
    return rng_;
#endif
}

// 函数说明：构造蒙特卡洛求解器并完成粒子、边界、热流统计与输出初始化。
MonteCarloSolver::MonteCarloSolver(const SimulationConfig& args, const SimulationDomain& geometry, const PhononMaterial& phonon)
    : args_(args), geometry_(&geometry), phonon_(&phonon) {
    using Clock = std::chrono::steady_clock;
    const auto t_init_begin = Clock::now();
    rng_seed_base_ = rng_();
    particle_count_ = std::max(1, static_cast<int>(std::llround(args_.particle_count)));
    time_step_ = std::max(1e-12, args_.time_step);
    push_eps_ = 1e-10 * std::max(time_step_, 1.0);
    particle_density_ = static_cast<double>(particle_count_) / std::max(geometry.volume(), 1e-12);

    std::cout << "MonteCarloSolver initialized: particle_count=" << particle_count_
              << ", time_step=" << time_step_
              << ", density=" << particle_density_ << '\n';
    std::cout << "Thermal conductivity estimation: "
              << (args_.compute_kappa ? "enabled" : "disabled")
              << '\n';
#ifdef _OPENMP
    openmp_thread_count_ = std::max(1, omp_get_max_threads());
    std::cout << "OpenMP enabled: max_threads=" << openmp_thread_count_ << '\n';
#else
    openmp_thread_count_ = 1;
    std::cout << "OpenMP enabled: no\n";
#endif

    if (!args_.output_folder.empty()) {
        std::filesystem::create_directories(args_.output_folder);
    }
    profile_timers_enabled_ = args_.profile_timers;
    if (profile_timers_enabled_) {
        std::cout << "Timestep profiling: enabled (simulation.profile_timers=true)\n";
    }

    initialize_particles(geometry, phonon);
    initialize_local_heat_source(geometry);
    write_convergence_header();
    update_heat_flux_and_conductivity(geometry);
    append_convergence_row();
    const auto t_init_end = Clock::now();
    const double init_sec = std::chrono::duration<double>(t_init_end - t_init_begin).count();
    std::cout << "[init] Initialization complete in " << std::fixed << std::setprecision(2)
              << init_sec << " s\n";
}

// 函数说明：初始化粒子主状态与碰撞缓存，建立时间推进的初始条件。
void MonteCarloSolver::initialize_particles(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    using Clock = std::chrono::steady_clock;
    auto step_begin = Clock::now();
    auto begin_step = [&](int idx, int total, const std::string& name) {
        step_begin = Clock::now();
        std::cout << "[init] " << idx << "/" << total << " " << name << "...\n";
    };
    auto end_step = [&]() {
        const auto step_end = Clock::now();
        const double sec = std::chrono::duration<double>(step_end - step_begin).count();
        std::cout << "        done (" << std::fixed << std::setprecision(2) << sec << " s)\n";
    };
    constexpr int total_steps = 8;

    const auto& mesh = geometry.mesh();
    begin_step(1, total_steps, "Sampling particle positions and assigning initial grid IDs");
    particle_positions_ = mesh.sample_volume_points(particle_count_);
    particle_grid_id_.assign(static_cast<size_t>(particle_count_), 0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        particle_grid_id_[i] = nearest_grid_index(geometry, particle_positions_[i]);
    }
    end_step();

    begin_step(2, total_steps, "Assigning phonon modes");
    initialize_particle_modes(phonon);
    end_step();
    begin_step(3, total_steps, "Initializing particle temperatures");
    initialize_particle_temperatures(geometry);
    end_step();
    begin_step(4, total_steps, "Initializing particle velocities");
    initialize_particle_velocities(phonon);
    end_step();
    begin_step(5, total_steps, "Building reservoir injection tables");
    initialize_reservoir_injection(geometry, phonon);
    end_step();
    begin_step(6, total_steps, "Precomputing rough-boundary scattering tables");
    initialize_rough_boundary_scattering(geometry, phonon);
    end_step();
    begin_step(7, total_steps, "Initializing particle state arrays and collision cache");
    particle_omega_.resize(static_cast<size_t>(particle_count_));
    particle_occupation_.resize(static_cast<size_t>(particle_count_));
    particle_energies_.resize(static_cast<size_t>(particle_count_));
    particle_alive_flags_.assign(static_cast<size_t>(particle_count_), static_cast<std::uint8_t>(1));
    const double tmin = particle_temperatures_.empty() ? 300.0 : *std::min_element(particle_temperatures_.begin(), particle_temperatures_.end());
    const int nsv = std::max(1, geometry.grid_count());
    grid_temperatures_.assign(static_cast<size_t>(nsv), tmin);
    grid_particle_counts_.assign(static_cast<size_t>(nsv), 0);
    grid_energy_density_.assign(static_cast<size_t>(nsv), 0.0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        particle_omega_[i] = phonon.mode_angular_frequency(particle_modes_[i]);
        particle_occupation_[i] = phonon.bose_occupation(particle_temperatures_[i], particle_modes_[i]);
        particle_energies_[i] = 0.0;
    }

    cached_collision_positions_.assign(static_cast<size_t>(particle_count_), {0.0, 0.0, 0.0});
    timesteps_to_collision_.assign(static_cast<size_t>(particle_count_), std::numeric_limits<double>::infinity());
    cached_collision_facets_.assign(static_cast<size_t>(particle_count_), -1);
    cached_collision_conditions_.assign(static_cast<size_t>(particle_count_), 'R');
    std::vector<int> all_idx(static_cast<size_t>(particle_count_));
    std::iota(all_idx.begin(), all_idx.end(), 0);
    update_collision_cache(geometry, all_idx);
    end_step();
    begin_step(8, total_steps, "Computing initial temperature field from particle energy");
    update_particle_temperatures(geometry, phonon);
    end_step();
}

// 函数说明：按材料活跃模态集合为每个粒子分配初始声子模态。
void MonteCarloSolver::initialize_particle_modes(const PhononMaterial& phonon) {
    particle_modes_.resize(static_cast<size_t>(particle_count_));
    for (int i = 0; i < particle_count_; ++i) {
        particle_modes_[i] = phonon.sample_active_mode(rng_);
    }
}

// 函数说明：依据边界温度与初始策略设置粒子温度场。
void MonteCarloSolver::initialize_particle_temperatures(const SimulationDomain& geometry) {
    particle_temperatures_.assign(static_cast<size_t>(particle_count_), 300.0);
    const auto& res_vals = geometry.reservoir_values();

    double tmin = 298.0;
    double tmax = 302.0;
    if (!res_vals.empty()) {
        tmin = *std::min_element(res_vals.begin(), res_vals.end());
        tmax = *std::max_element(res_vals.begin(), res_vals.end());
        if (std::isnan(tmin) || std::isnan(tmax)) {
            tmin = 298.0;
            tmax = 302.0;
        }
    }
    if (tmin > tmax) {
        std::swap(tmin, tmax);
    }
    const double tmean = 0.5 * (tmin + tmax);
    std::string key = "cold";
    if (!args_.initial_temperature.empty()) {
        key = args_.initial_temperature.front();
    }

    if (key == "cold") {
        std::fill(particle_temperatures_.begin(), particle_temperatures_.end(), tmin);
    } else if (key == "mean") {
        std::fill(particle_temperatures_.begin(), particle_temperatures_.end(), tmean);
    } else {
        throw std::runtime_error("initial_temperature only supports: cold, mean");
    }
}

// 函数说明：采样随机单位方向，用于漫反射与方向随机化过程。
MonteCarloSolver::Vec3 MonteCarloSolver::random_unit_vector() {
    auto& rng = thread_rng();
    std::normal_distribution<double> N(0.0, 1.0);
    Vec3 v {N(rng), N(rng), N(rng)};
    const double n = norm(v);
    if (n <= 1e-12) {
        return {1.0, 0.0, 0.0};
    }
    return mul(v, 1.0 / n);
}

// 函数说明：根据粒子模态查询群速度并写入速度场。
void MonteCarloSolver::initialize_particle_velocities(const PhononMaterial& phonon) {
    particle_velocities_.resize(static_cast<size_t>(particle_count_));
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        particle_velocities_[i] = phonon.mode_group_velocity(particle_modes_[i]);
    }
}

// 函数说明：预计算热库注入概率与计数器，驱动边界粒子注入机制。
void MonteCarloSolver::initialize_reservoir_injection(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    reservoir_facets_ = geometry.reservoir_facets();
    reservoir_temperatures_.clear();
    reservoir_temperatures_.reserve(reservoir_facets_.size());
    for (int facet : reservoir_facets_) {
        reservoir_temperatures_.push_back(geometry.reservoir_value_for_facet(facet, 300.0));
    }
    reservoir_count_ = static_cast<int>(reservoir_facets_.size());
    reservoir_modes_ = phonon.active_mode_list();
    reservoir_entry_probability_.assign(static_cast<size_t>(reservoir_count_), std::vector<double>(reservoir_modes_.size(), 0.0));
    reservoir_emit_counter_.assign(static_cast<size_t>(reservoir_count_), std::vector<double>(reservoir_modes_.size(), 0.0));
    reservoir_areas_.assign(static_cast<size_t>(reservoir_count_), 0.0);
    reservoir_normals_.assign(static_cast<size_t>(reservoir_count_), {0.0, 0.0, 1.0});

    auto& rng = thread_rng();
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const auto& mesh = geometry.mesh();
    const auto& areas = mesh.facet_areas();
    const auto& normals = mesh.facet_normals();
    const double n_active = std::max(1.0, static_cast<double>(phonon.active_mode_count()));

    for (int r = 0; r < reservoir_count_; ++r) {
        const int facet = reservoir_facets_[static_cast<size_t>(r)];
        if (facet < 0 || facet >= static_cast<int>(areas.size()) || facet >= static_cast<int>(normals.size())) {
            continue;
        }
        reservoir_areas_[static_cast<size_t>(r)] = std::max(1e-12, areas[static_cast<size_t>(facet)]);
        reservoir_normals_[static_cast<size_t>(r)] = normals[static_cast<size_t>(facet)];
        const double bound_thickness = n_active / std::max(1e-18, particle_density_ * reservoir_areas_[static_cast<size_t>(r)]);
        for (size_t m = 0; m < reservoir_modes_.size(); ++m) {
            const Vec3 gv = phonon.mode_group_velocity(reservoir_modes_[m]);
            const double gv_parallel = -dot(reservoir_normals_[static_cast<size_t>(r)], gv);
            const double p = std::max(0.0, gv_parallel * time_step_ / std::max(1e-18, bound_thickness));
            reservoir_entry_probability_[static_cast<size_t>(r)][m] = p;
            reservoir_emit_counter_[static_cast<size_t>(r)][m] = U(rng);
        }
    }
}

// 函数说明：构建粗糙边界散射查找表，包含镜面率与漫反射候选映射。
void MonteCarloSolver::initialize_rough_boundary_scattering(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    const auto& mesh = geometry.mesh();
    const int nfacets = mesh.facet_count();
    facet_to_rough_data_.assign(static_cast<size_t>(nfacets), -1);
    rough_boundary_data_.clear();

    const int na = phonon.active_mode_count();
    if (na <= 0) {
        return;
    }
    const auto& active = phonon.active_mode_list();
    std::vector<Vec3> v_active(static_cast<size_t>(na), {0.0, 0.0, 0.0});
    std::vector<double> vnorm_active(static_cast<size_t>(na), 0.0);
    std::vector<double> omega_active(static_cast<size_t>(na), 0.0);
    std::vector<double> knorm_active(static_cast<size_t>(na), 0.0);
    std::vector<double> domega_active(static_cast<size_t>(na), 0.0);
    for (int ai = 0; ai < na; ++ai) {
        v_active[static_cast<size_t>(ai)] = phonon.mode_group_velocity(active[static_cast<size_t>(ai)]);
        vnorm_active[static_cast<size_t>(ai)] = norm(v_active[static_cast<size_t>(ai)]);
        omega_active[static_cast<size_t>(ai)] = phonon.mode_angular_frequency(active[static_cast<size_t>(ai)]);
        knorm_active[static_cast<size_t>(ai)] = phonon.mode_wavevector_norm(active[static_cast<size_t>(ai)]);
        domega_active[static_cast<size_t>(ai)] = phonon.mode_frequency_window(active[static_cast<size_t>(ai)]);
    }

    int rough_total = 0;
    for (int facet = 0; facet < nfacets; ++facet) {
        if (geometry.is_rough_facet(facet)) {
            ++rough_total;
        }
    }
    if (rough_total <= 0) {
        std::cout << "[init] Rough-boundary preprocessing skipped (0 rough facets).\n";
        return;
    }
    int rough_done = 0;
    const int report_stride = std::max(1, rough_total / 10);

    for (int facet = 0; facet < nfacets; ++facet) {
        if (!geometry.is_rough_facet(facet)) {
            continue;
        }
        ++rough_done;
        if (rough_done == 1 || rough_done == rough_total || (rough_done % report_stride) == 0) {
            std::cout << "[init]   rough facet " << rough_done << "/" << rough_total
                      << " (mesh facet id=" << facet << ")\n";
        }
        RoughFacetData rd;
        rd.facet = facet;
        rd.specularity.assign(static_cast<size_t>(na), 0.0);
        rd.spec_match_active.assign(static_cast<size_t>(na), -1);
        rd.diffuse_creation_rate.assign(static_cast<size_t>(na), 0.0);
        rd.diffuse_creation_prob.assign(static_cast<size_t>(na), 0.0);
        rd.diffuse_begin.assign(static_cast<size_t>(na), -1);
        rd.diffuse_end.assign(static_cast<size_t>(na), -1);

        const double eta = std::max(0.0, geometry.roughness_for_facet(facet, 0.0));
        const Vec3 n_out = mesh.facet_normals()[static_cast<size_t>(facet)];
        const Vec3 n_in {-n_out[0], -n_out[1], -n_out[2]};
        std::vector<int> incoming;
        std::vector<int> outgoing;
        std::vector<double> incoming_destruction(static_cast<size_t>(na), 0.0);
        std::vector<double> outgoing_creation(static_cast<size_t>(na), 0.0);
        incoming.reserve(static_cast<size_t>(na));
        outgoing.reserve(static_cast<size_t>(na));

        for (int ai = 0; ai < na; ++ai) {
            const Vec3 v = v_active[static_cast<size_t>(ai)];
            const double vn = dot(v, n_in);
            if (vn < 0.0) {
                incoming.push_back(ai);
                incoming_destruction[static_cast<size_t>(ai)] = -vn;
            } else if (vn > 0.0) {
                outgoing.push_back(ai);
                outgoing_creation[static_cast<size_t>(ai)] = vn;
            }
            const double vnorm = std::max(1e-12, vnorm_active[static_cast<size_t>(ai)]);
            const double incidence_cos = std::abs(vn) / vnorm;
            const double x = 2.0 * eta * incidence_cos * knorm_active[static_cast<size_t>(ai)];
            double p = std::exp(-(x * x));
            if (!std::isfinite(p)) {
                p = 0.0;
            }
            // Specularity is only meaningful for incoming modes on this facet.
            rd.specularity[static_cast<size_t>(ai)] = (vn < 0.0) ? std::clamp(p, 0.0, 1.0) : 0.0;
        }
        rd.outgoing_active = outgoing;
        rd.outgoing_sorted_active = outgoing;
        std::sort(
            rd.outgoing_sorted_active.begin(),
            rd.outgoing_sorted_active.end(),
            [&omega_active](int a, int b) { return omega_active[static_cast<size_t>(a)] < omega_active[static_cast<size_t>(b)]; });
        rd.outgoing_sorted_omega.resize(rd.outgoing_sorted_active.size());
        for (size_t k = 0; k < rd.outgoing_sorted_active.size(); ++k) {
            rd.outgoing_sorted_omega[k] = omega_active[static_cast<size_t>(rd.outgoing_sorted_active[k])];
        }

        for (int ai : incoming) {
            if (rd.outgoing_sorted_active.empty()) {
                continue;
            }
            const double w_in = omega_active[static_cast<size_t>(ai)];
            const double tol = std::max(1e-12, domega_active[static_cast<size_t>(ai)]);
            const double w_lo = w_in - tol;
            const double w_hi = w_in + tol;
            auto it_lo = std::lower_bound(rd.outgoing_sorted_omega.begin(), rd.outgoing_sorted_omega.end(), w_lo);
            auto it_hi = std::upper_bound(rd.outgoing_sorted_omega.begin(), rd.outgoing_sorted_omega.end(), w_hi);
            const int ib = static_cast<int>(std::distance(rd.outgoing_sorted_omega.begin(), it_lo));
            const int ie = static_cast<int>(std::distance(rd.outgoing_sorted_omega.begin(), it_hi));
            if (ie <= ib) {
                continue;
            }
            rd.diffuse_begin[static_cast<size_t>(ai)] = ib;
            rd.diffuse_end[static_cast<size_t>(ai)] = ie;

            const Vec3 vin = v_active[static_cast<size_t>(ai)];
            const double vin_dot_n = dot(vin, n_in);
            const Vec3 vtry = sub(vin, mul(n_in, 2.0 * vin_dot_n));

            int best_ao = -1;
            double best_metric = std::numeric_limits<double>::infinity();
            for (int pos = ib; pos < ie; ++pos) {
                const int ao = rd.outgoing_sorted_active[static_cast<size_t>(pos)];
                const Vec3 vout = v_active[static_cast<size_t>(ao)];
                const Vec3 dv = sub(vtry, vout);
                const double refn = std::max({1e-12, norm(vtry), norm(vout)});
                const double rel = norm(dv) / refn;
                const double metric = rel;
                if (metric < best_metric) {
                    best_metric = metric;
                    best_ao = ao;
                }
            }
            if (best_ao >= 0) {
                rd.spec_match_active[static_cast<size_t>(ai)] = best_ao;
            }
        }

        // Population.py style: only truly matched incoming modes can be specular.
        for (int ai : incoming) {
            if (rd.spec_match_active[static_cast<size_t>(ai)] < 0) {
                rd.specularity[static_cast<size_t>(ai)] = 0.0;
            }
        }

        // Population.py style diffuse creation rate:
        // creation_rate(out) = C_total(out) - sum(specular_D(in -> out)).
        for (int ao : outgoing) {
            rd.diffuse_creation_rate[static_cast<size_t>(ao)] = outgoing_creation[static_cast<size_t>(ao)];
        }
        for (int ai : incoming) {
            const int ao = rd.spec_match_active[static_cast<size_t>(ai)];
            if (ao < 0 || ao >= na) {
                continue;
            }
            const double specular_d =
                incoming_destruction[static_cast<size_t>(ai)] * rd.specularity[static_cast<size_t>(ai)];
            rd.diffuse_creation_rate[static_cast<size_t>(ao)] -= specular_d;
        }

        double rate_sum = 0.0;
        for (int ao : outgoing) {
            double& r = rd.diffuse_creation_rate[static_cast<size_t>(ao)];
            if (!std::isfinite(r) || r < 0.0) {
                r = 0.0;
            }
            rate_sum += r;
        }
        if (rate_sum > 0.0) {
            double cdf = 0.0;
            for (int ao : outgoing) {
                const double r = rd.diffuse_creation_rate[static_cast<size_t>(ao)];
                if (r <= 0.0) {
                    continue;
                }
                cdf += r;
                rd.diffuse_roulette_active.push_back(ao);
                rd.diffuse_roulette_cdf.push_back(cdf);
                rd.diffuse_creation_prob[static_cast<size_t>(ao)] = r / rate_sum;
            }
        }

        const int rid = static_cast<int>(rough_boundary_data_.size());
        rough_boundary_data_.push_back(std::move(rd));
        facet_to_rough_data_[static_cast<size_t>(facet)] = rid;
    }
}

// 函数说明：导出粗糙面模态映射，便于核查镜面/漫反射模式选择是否按预计算执行。
void MonteCarloSolver::write_rough_boundary_mode_map(const SimulationDomain& geometry, const PhononMaterial& phonon) const {
    if (args_.output_folder.empty()) {
        return;
    }
    namespace fs = std::filesystem;
    fs::create_directories(args_.output_folder);
    std::ofstream out(fs::path(args_.output_folder) / "rough_boundary_mode_map.csv", std::ios::trunc);
    if (!out) {
        return;
    }

    out << "rough_index,facet,in_active_idx,in_q,in_branch,p_spec,spec_match_active,spec_match_q,spec_match_branch,"
           "diffuse_begin,diffuse_end,diffuse_candidate_count,outgoing_count,is_incoming,creation_rate,creation_prob,diffuse_roulette_size\n";

    const int na = phonon.active_mode_count();
    for (int rid = 0; rid < static_cast<int>(rough_boundary_data_.size()); ++rid) {
        const auto& rd = rough_boundary_data_[static_cast<size_t>(rid)];
        for (int ai = 0; ai < na; ++ai) {
            const auto in_mode = phonon.active_mode_at(ai);
            const double p_spec = (ai < static_cast<int>(rd.specularity.size())) ? rd.specularity[static_cast<size_t>(ai)] : 0.0;
            const int spec_ai = (ai < static_cast<int>(rd.spec_match_active.size())) ? rd.spec_match_active[static_cast<size_t>(ai)] : -1;
            const auto spec_mode = (spec_ai >= 0 && spec_ai < na) ? phonon.active_mode_at(spec_ai) : std::array<int, 2>{-1, -1};
            const int ib = (ai < static_cast<int>(rd.diffuse_begin.size())) ? rd.diffuse_begin[static_cast<size_t>(ai)] : -1;
            const int ie = (ai < static_cast<int>(rd.diffuse_end.size())) ? rd.diffuse_end[static_cast<size_t>(ai)] : -1;
            const int n_candidate = (ib >= 0 && ie >= ib) ? (ie - ib) : 0;
            const int is_incoming = (n_candidate > 0) ? 1 : 0;
            const double creation_rate =
                (ai < static_cast<int>(rd.diffuse_creation_rate.size())) ? rd.diffuse_creation_rate[static_cast<size_t>(ai)] : 0.0;
            const double creation_prob =
                (ai < static_cast<int>(rd.diffuse_creation_prob.size())) ? rd.diffuse_creation_prob[static_cast<size_t>(ai)] : 0.0;

            out << rid << ","
                << rd.facet << ","
                << ai << ","
                << in_mode[0] << ","
                << in_mode[1] << ","
                << std::setprecision(17) << p_spec << ","
                << spec_ai << ","
                << spec_mode[0] << ","
                << spec_mode[1] << ","
                << ib << ","
                << ie << ","
                << n_candidate << ","
                << rd.outgoing_active.size() << ","
                << is_incoming << ","
                << creation_rate << ","
                << creation_prob << ","
                << rd.diffuse_roulette_active.size() << "\n";
        }
    }

    std::cout << "Rough-boundary mode map written: "
              << (fs::path(args_.output_folder) / "rough_boundary_mode_map.csv").string() << '\n';
    (void) geometry;
}

// 函数说明：解析局部热源区域并生成网格掩码与功率参数。
void MonteCarloSolver::initialize_local_heat_source(const SimulationDomain& geometry) {
    local_heat_source_enabled_ = false;
    local_heat_source_grid_mask_.clear();
    local_heat_source_power_density_wm3_ = args_.heat_source_power_density;

    if (!args_.heat_source_enabled) {
        return;
    }
    if (args_.heat_source_min.size() != 3 || args_.heat_source_max.size() != 3) {
        std::cerr << "Warning: heat source enabled but min/max are not 3D vectors. Ignoring local heat source.\n";
        return;
    }
    if (std::abs(local_heat_source_power_density_wm3_) <= 0.0) {
        std::cerr << "Warning: heat source enabled but power_density is zero. Ignoring local heat source.\n";
        return;
    }

    local_heat_source_min_ = {
        args_.heat_source_min[0], args_.heat_source_min[1], args_.heat_source_min[2]
    };
    local_heat_source_max_ = {
        args_.heat_source_max[0], args_.heat_source_max[1], args_.heat_source_max[2]
    };

    // Heat source coordinates are always interpreted as relative [0, 1] box fractions.
    const auto& bmin = geometry.bounds_min();
    const auto& bmax = geometry.bounds_max();
    const Vec3 ext {bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]};
    for (int k = 0; k < 3; ++k) {
        local_heat_source_min_[k] = bmin[k] + local_heat_source_min_[k] * ext[k];
        local_heat_source_max_[k] = bmin[k] + local_heat_source_max_[k] * ext[k];
    }

    for (int k = 0; k < 3; ++k) {
        if (local_heat_source_min_[k] > local_heat_source_max_[k]) {
            std::swap(local_heat_source_min_[k], local_heat_source_max_[k]);
        }
    }

    const auto& centers = geometry.grid_centers();
    local_heat_source_grid_mask_.assign(centers.size(), static_cast<std::uint8_t>(0));
    int selected = 0;
    for (size_t i = 0; i < centers.size(); ++i) {
        const Vec3 c = centers[i];
        const bool inside =
            c[0] >= local_heat_source_min_[0] && c[0] <= local_heat_source_max_[0] &&
            c[1] >= local_heat_source_min_[1] && c[1] <= local_heat_source_max_[1] &&
            c[2] >= local_heat_source_min_[2] && c[2] <= local_heat_source_max_[2];
        if (inside) {
            local_heat_source_grid_mask_[i] = static_cast<std::uint8_t>(1);
            ++selected;
        }
    }

    if (selected == 0) {
        std::cerr << "Warning: heat source region does not include any grid center. Ignoring local heat source.\n";
        return;
    }

    local_heat_source_enabled_ = true;
    std::cout << "Local heat source enabled: selected_grids=" << selected
              << ", power_density=" << local_heat_source_power_density_wm3_ << " W/m^3\n";
}

// 函数说明：在每步能量更新后向热源区域叠加体热源能量。
void MonteCarloSolver::apply_local_heat_source() {
    if (!local_heat_source_enabled_) {
        return;
    }
    const double delta_e = local_heat_source_power_density_wm3_ * wm3_to_evpsa3_ * time_step_;
    if (!std::isfinite(delta_e) || std::abs(delta_e) <= 0.0) {
        return;
    }
    const size_t n = std::min(grid_energy_density_.size(), local_heat_source_grid_mask_.size());
    for (size_t i = 0; i < n; ++i) {
        if (local_heat_source_grid_mask_[i] != 0) {
            grid_energy_density_[i] += delta_e;
        }
    }
}

// 函数说明：在粗糙边界条件下采样漫反射后的活跃模态索引。
// source: 1=creation-roulette, 2=outgoing-pool fallback, 3=global-random fallback.
int MonteCarloSolver::sample_diffuse_active_mode(int rough_idx, int in_ai, int* source) const {
    auto& rng = thread_rng();
    const int na = (phonon_ != nullptr) ? phonon_->active_mode_count() : 0;
    (void) in_ai;
    if (na <= 0) {
        if (source != nullptr) {
            *source = 3;
        }
        return 0;
    }
    if (rough_idx < 0 || rough_idx >= static_cast<int>(rough_boundary_data_.size())) {
        std::uniform_int_distribution<int> U(0, na - 1);
        if (source != nullptr) {
            *source = 3;
        }
        return U(rng);
    }
    const auto& rd = rough_boundary_data_[static_cast<size_t>(rough_idx)];

    if (!rd.diffuse_roulette_active.empty() &&
        !rd.diffuse_roulette_cdf.empty() &&
        rd.diffuse_roulette_active.size() == rd.diffuse_roulette_cdf.size() &&
        rd.diffuse_roulette_cdf.back() > 0.0) {
        std::uniform_real_distribution<double> U(0.0, rd.diffuse_roulette_cdf.back());
        const double r = U(rng);
        auto it = std::lower_bound(rd.diffuse_roulette_cdf.begin(), rd.diffuse_roulette_cdf.end(), r);
        size_t pos = static_cast<size_t>(std::distance(rd.diffuse_roulette_cdf.begin(), it));
        if (pos >= rd.diffuse_roulette_active.size()) {
            pos = rd.diffuse_roulette_active.size() - 1;
        }
        if (source != nullptr) {
            *source = 1;
        }
        return rd.diffuse_roulette_active[pos];
    }

    if (!rd.outgoing_active.empty()) {
        std::uniform_int_distribution<int> U(0, static_cast<int>(rd.outgoing_active.size()) - 1);
        if (source != nullptr) {
            *source = 2;
        }
        return rd.outgoing_active[static_cast<size_t>(U(rng))];
    }
    std::uniform_int_distribution<int> U(0, na - 1);
    if (source != nullptr) {
        *source = 3;
    }
    return U(rng);
}

// 函数说明：按镜面/漫反射概率选择碰撞后的模态与占据数。
std::array<int, 2> MonteCarloSolver::select_reflected_mode(
    const SimulationDomain& geometry,
    const PhononMaterial& phonon,
    int rough_idx,
    const std::array<int, 2>& in_mode,
    const Vec3& collision_pos,
    double& out_occupation,
    double in_occupation) const {
    const int na = phonon.active_mode_count();
    rough_events_total_.fetch_add(1, std::memory_order_relaxed);
    if (rough_idx < 0 || rough_idx >= static_cast<int>(rough_boundary_data_.size()) || na <= 0) {
        rough_fallback_missing_rough_data_.fetch_add(1, std::memory_order_relaxed);
        out_occupation = in_occupation;
        return in_mode;
    }
    const auto& rd = rough_boundary_data_[static_cast<size_t>(rough_idx)];
    int in_ai = phonon.active_index_for_mode(in_mode);
    if (in_ai < 0 || in_ai >= na) {
        in_ai = -1;
    }

    std::uniform_real_distribution<double> U01(0.0, 1.0);
    auto& rng = thread_rng();
    const double p_spec = (in_ai >= 0 && in_ai < static_cast<int>(rd.specularity.size()))
        ? std::clamp(rd.specularity[static_cast<size_t>(in_ai)], 0.0, 1.0)
        : 0.0;
    if (in_ai >= 0 && U01(rng) <= p_spec) {
        const int out_ai = rd.spec_match_active[static_cast<size_t>(in_ai)];
        if (out_ai >= 0 && out_ai < na) {
            rough_specular_selected_.fetch_add(1, std::memory_order_relaxed);
            out_occupation = in_occupation;
            return phonon.active_mode_at(out_ai);
        }
        rough_fallback_missing_spec_match_.fetch_add(1, std::memory_order_relaxed);
    }

    int source = 3;
    const int out_ai = sample_diffuse_active_mode(rough_idx, in_ai, &source);
    rough_diffuse_selected_.fetch_add(1, std::memory_order_relaxed);
    if (source == 2) {
        rough_fallback_outgoing_pool_.fetch_add(1, std::memory_order_relaxed);
    } else if (source == 3) {
        rough_fallback_global_random_.fetch_add(1, std::memory_order_relaxed);
    }
    const std::array<int, 2> out_mode = phonon.active_mode_at(out_ai);
    const int nsv = std::max(1, geometry.grid_count());
    const int sv = std::clamp(nearest_grid_index(geometry, collision_pos), 0, nsv - 1);
    const double Tdiff = (sv >= 0 && sv < static_cast<int>(grid_temperatures_.size()))
        ? grid_temperatures_[static_cast<size_t>(sv)] : 300.0;
    out_occupation = phonon.bose_occupation(Tdiff, out_mode);
    return out_mode;
}

// 函数说明：将粒子位置映射到最近控制体网格索引。
int MonteCarloSolver::nearest_grid_index(const SimulationDomain& geometry, const Vec3& p) const {
    const auto& centers = geometry.grid_centers();
    if (centers.empty()) {
        return 0;
    }
    int best = 0;
    double best_d2 = std::numeric_limits<double>::max();
    for (int i = 0; i < static_cast<int>(centers.size()); ++i) {
        const auto d = sub(p, centers[i]);
        const double d2 = dot(d, d);
        if (d2 < best_d2) {
            best_d2 = d2;
            best = i;
        }
    }
    return best;
}

// 函数说明：批量追踪粒子到下一次边界碰撞点并缓存结果。
void MonteCarloSolver::update_collision_cache(const SimulationDomain& geometry, const std::vector<int>& indices) {
    const int nidx = static_cast<int>(indices.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int k = 0; k < nidx; ++k) {
        update_collision_cache_single(geometry, indices[static_cast<size_t>(k)]);
    }
}

// 函数说明：为单个粒子更新下一次边界碰撞缓存。
void MonteCarloSolver::update_collision_cache_single(const SimulationDomain& geometry, int i) {
    const auto& mesh = geometry.mesh();
    if (i < 0 || i >= particle_count_) {
        return;
    }
    auto [cp, t_hit, fct] = mesh.trace_boundary_intersection(particle_positions_[i], particle_velocities_[i]);
    if (fct < 0 || std::isinf(t_hit)) {
        // Robust fallback: reverse direction and retry once.
        particle_velocities_[i] = mul(particle_velocities_[i], -1.0);
        std::tie(cp, t_hit, fct) = mesh.trace_boundary_intersection(particle_positions_[i], particle_velocities_[i]);
    }
    cached_collision_positions_[i] = cp;
    cached_collision_facets_[i] = fct;
    cached_collision_conditions_[i] = geometry.facet_boundary_condition(fct);
    timesteps_to_collision_[i] = std::isinf(t_hit) ? std::numeric_limits<double>::infinity() : (t_hit / time_step_);
}

// 函数说明：处理边界事件（透射、吸收、周期、粗糙/镜面反射）并更新粒子态。
void MonteCarloSolver::process_boundary_collision(const SimulationDomain& geometry, int i) {
    const int facet = cached_collision_facets_[i];
    const char cond = cached_collision_conditions_[i];

    particle_positions_[i] = cached_collision_positions_[i];
    if (facet < 0 || facet >= static_cast<int>(geometry.mesh().facet_normals().size())) {
        return;
    }

    const Vec3 n = geometry.mesh().facet_normals()[facet];
    const double vn = dot(particle_velocities_[i], n);

    if (cond == 'T') {
        if (i >= 0 && i < static_cast<int>(particle_alive_flags_.size())) {
            particle_alive_flags_[static_cast<size_t>(i)] = static_cast<std::uint8_t>(0);
        }
        return;
    } else if (cond == 'P' && geometry.has_periodic_pair(facet)) {
        // Teleport through periodic interface with facet-pair translation.
        const Vec3 shift = geometry.periodic_shift_for_facet(facet);
        particle_positions_[i] = add(particle_positions_[i], shift);
    } else if (cond == 'R' && phonon_ != nullptr) {
        const int rough_idx = (facet >= 0 && facet < static_cast<int>(facet_to_rough_data_.size()))
            ? facet_to_rough_data_[static_cast<size_t>(facet)] : -1;
        if (rough_idx >= 0) {
            const double in_occ = particle_occupation_[i];
            double out_occ = in_occ;
            const std::array<int, 2> out_mode =
                select_reflected_mode(geometry, *phonon_, rough_idx, particle_modes_[i], cached_collision_positions_[i], out_occ, in_occ);
            particle_modes_[i] = out_mode;
            particle_velocities_[i] = phonon_->mode_group_velocity(particle_modes_[i]);
            if (dot(particle_velocities_[i], n) > 0.0) {
                particle_velocities_[i] = sub(particle_velocities_[i], mul(n, 2.0 * dot(particle_velocities_[i], n)));
            }
            particle_omega_[i] = phonon_->mode_angular_frequency(particle_modes_[i]);
            particle_occupation_[i] = out_occ;
            particle_energies_[i] = 0.0;
        } else {
            rough_events_total_.fetch_add(1, std::memory_order_relaxed);
            rough_fallback_missing_rough_data_.fetch_add(1, std::memory_order_relaxed);
            const double p_spec = compute_roughness_specularity(geometry, *phonon_, i, facet);
            std::uniform_real_distribution<double> U(0.0, 1.0);
            auto& rng = thread_rng();
            if (U(rng) <= p_spec) {
                rough_specular_selected_.fetch_add(1, std::memory_order_relaxed);
                particle_velocities_[i] = sub(particle_velocities_[i], mul(n, 2.0 * vn));
            } else {
                rough_diffuse_selected_.fetch_add(1, std::memory_order_relaxed);
                rough_fallback_global_random_.fetch_add(1, std::memory_order_relaxed);
                particle_modes_[i] = phonon_->sample_active_mode(rng);
                const double speed = std::max(1e-9, norm(phonon_->mode_group_velocity(particle_modes_[i])));
                Vec3 dir = random_unit_vector();
                if (dot(dir, n) > 0.0) {
                    dir = mul(dir, -1.0);
                }
                particle_velocities_[i] = mul(dir, speed);
                particle_omega_[i] = phonon_->mode_angular_frequency(particle_modes_[i]);
                particle_occupation_[i] = phonon_->bose_occupation(particle_temperatures_[i], particle_modes_[i]);
                particle_energies_[i] = 0.0;
            }
        }
    } else {
        // Unmatched periodic and generic boundaries: specular reflection.
        particle_velocities_[i] = sub(particle_velocities_[i], mul(n, 2.0 * vn));
    }

    particle_positions_[i] = add(particle_positions_[i], mul(particle_velocities_[i], push_eps_ / std::max(norm(particle_velocities_[i]), 1e-12)));
    particle_grid_id_[i] = nearest_grid_index(geometry, particle_positions_[i]);
}

// 函数说明：根据粗糙度与入射条件计算镜面反射概率。
double MonteCarloSolver::compute_roughness_specularity(const SimulationDomain& geometry, const PhononMaterial& phonon, int i, int facet) const {
    const double eta = std::max(0.0, geometry.roughness_for_facet(facet, 0.0));
    const double speed = std::max(1e-12, norm(particle_velocities_[i]));
    const Vec3 n = geometry.mesh().facet_normals()[facet];
    const double incidence_cos = std::min(1.0, std::max(0.0, std::abs(dot(particle_velocities_[i], n)) / speed));
    const double k_norm = std::max(1e-12, phonon.mode_angular_frequency(particle_modes_[i]) / speed);
    const double x = 2.0 * eta * incidence_cos * k_norm;
    const double p = std::exp(-(x * x));
    return std::clamp(p, 0.0, 1.0);
}

// 函数说明：移除已吸收粒子并压缩所有粒子状态数组。
void MonteCarloSolver::remove_absorbed_particles() {
    if (particle_alive_flags_.empty()) {
        return;
    }
    size_t alive_count = 0;
    for (std::uint8_t a : particle_alive_flags_) {
        alive_count += (a != 0) ? 1u : 0u;
    }
    if (alive_count == particle_alive_flags_.size()) {
        return;
    }

    std::vector<size_t> keep;
    keep.reserve(alive_count);
    for (size_t i = 0; i < particle_alive_flags_.size(); ++i) {
        if (particle_alive_flags_[i] != 0) {
            keep.push_back(i);
        }
    }

    auto remap_vec3 = [&keep](std::vector<Vec3>& v) {
        std::vector<Vec3> out;
        out.reserve(keep.size());
        for (size_t k : keep) {
            out.push_back(v[k]);
        }
        v.swap(out);
    };
    auto remap_vecd = [&keep](std::vector<double>& v) {
        std::vector<double> out;
        out.reserve(keep.size());
        for (size_t k : keep) {
            out.push_back(v[k]);
        }
        v.swap(out);
    };
    auto remap_veci = [&keep](std::vector<int>& v) {
        std::vector<int> out;
        out.reserve(keep.size());
        for (size_t k : keep) {
            out.push_back(v[k]);
        }
        v.swap(out);
    };
    auto remap_vecc = [&keep](std::vector<char>& v) {
        std::vector<char> out;
        out.reserve(keep.size());
        for (size_t k : keep) {
            out.push_back(v[k]);
        }
        v.swap(out);
    };
    auto remap_modes = [&keep](std::vector<std::array<int, 2>>& v) {
        std::vector<std::array<int, 2>> out;
        out.reserve(keep.size());
        for (size_t k : keep) {
            out.push_back(v[k]);
        }
        v.swap(out);
    };

    remap_modes(particle_modes_);
    remap_vec3(particle_positions_);
    remap_vec3(particle_velocities_);
    remap_vec3(cached_collision_positions_);
    remap_vecd(timesteps_to_collision_);
    remap_vecd(particle_temperatures_);
    remap_vecd(particle_omega_);
    remap_vecd(particle_occupation_);
    remap_vecd(particle_energies_);
    remap_veci(particle_grid_id_);
    remap_veci(cached_collision_facets_);
    remap_vecc(cached_collision_conditions_);
    particle_alive_flags_.assign(alive_count, static_cast<std::uint8_t>(1));
    particle_count_ = static_cast<int>(alive_count);
}

// 函数说明：将数值误差导致越界的粒子重采样回几何体内部，并重建碰撞缓存。
void MonteCarloSolver::recover_escaped_particles(const SimulationDomain& geometry) {
    if (particle_count_ <= 0 || particle_positions_.empty()) {
        return;
    }
    const auto& bmin = geometry.bounds_min();
    const auto& bmax = geometry.bounds_max();
    const Vec3 ext {bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]};
    const double scale = std::max({1.0, std::abs(ext[0]), std::abs(ext[1]), std::abs(ext[2])});
    const double tol = 1e-10 * scale;

    std::vector<int> escaped_idx;
    escaped_idx.reserve(64);
    for (int i = 0; i < particle_count_; ++i) {
        const Vec3& p = particle_positions_[static_cast<size_t>(i)];
        const bool out =
            (p[0] < bmin[0] - tol || p[0] > bmax[0] + tol) ||
            (p[1] < bmin[1] - tol || p[1] > bmax[1] + tol) ||
            (p[2] < bmin[2] - tol || p[2] > bmax[2] + tol);
        if (out) {
            escaped_idx.push_back(i);
        }
    }
    if (escaped_idx.empty()) {
        return;
    }

    const auto& mesh = geometry.mesh();
    std::vector<Vec3> respawn = mesh.sample_volume_points(static_cast<int>(escaped_idx.size()));
    for (size_t k = 0; k < escaped_idx.size(); ++k) {
        const int idx = escaped_idx[k];
        particle_positions_[static_cast<size_t>(idx)] = respawn[k];
        particle_grid_id_[static_cast<size_t>(idx)] = nearest_grid_index(geometry, respawn[k]);
    }
    update_collision_cache(geometry, escaped_idx);

    escaped_recovery_events_ += 1;
    escaped_recovered_particles_ += static_cast<long long>(escaped_idx.size());
}

// 函数说明：从热库边界按统计规则注入新粒子并返回注入时刻偏移。
std::vector<std::pair<int, double>> MonteCarloSolver::inject_particles_from_reservoirs(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    std::vector<std::pair<int, double>> inserted;
    if (reservoir_count_ <= 0 || reservoir_modes_.empty()) {
        return inserted;
    }

    auto& rng = thread_rng();
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const auto& mesh = geometry.mesh();

    for (int r = 0; r < reservoir_count_; ++r) {
        const int facet = reservoir_facets_[static_cast<size_t>(r)];
        std::vector<size_t> mode_idx;
        std::vector<double> dt_in;
        mode_idx.reserve(128);
        dt_in.reserve(128);
        auto& row_prob = reservoir_entry_probability_[static_cast<size_t>(r)];
        auto& row_counter = reservoir_emit_counter_[static_cast<size_t>(r)];

        for (size_t m = 0; m < row_prob.size(); ++m) {
            const double p = std::max(0.0, row_prob[m]);
            if (p <= 0.0) {
                continue;
            }
            const int fixed = static_cast<int>(std::floor(p));
            row_counter[m] += p - static_cast<double>(fixed);
            int n_emit = fixed;
            while (row_counter[m] >= 1.0) {
                row_counter[m] -= 1.0;
                ++n_emit;
            }
            for (int c = 0; c < n_emit; ++c) {
                mode_idx.push_back(m);
                dt_in.push_back(time_step_ * U01(rng));
            }
        }

        if (mode_idx.empty()) {
            continue;
        }
        std::vector<Vec3> surf = mesh.sample_surface_points(static_cast<int>(mode_idx.size()), std::vector<int>{facet});
        const Vec3 n = reservoir_normals_[static_cast<size_t>(r)];
        const double Tres = reservoir_temperatures_[static_cast<size_t>(r)];

        for (size_t k = 0; k < mode_idx.size(); ++k) {
            const std::array<int, 2> mode = reservoir_modes_[mode_idx[k]];
            Vec3 gv = phonon.mode_group_velocity(mode);
            if (dot(gv, n) > 0.0) {
                gv = mul(gv, -1.0);
            }
            const double speed = std::max(1e-12, norm(gv));
            Vec3 pos = add(surf[k], mul(gv, push_eps_ / speed));

            const int new_idx = particle_count_;
            particle_modes_.push_back(mode);
            particle_positions_.push_back(pos);
            particle_velocities_.push_back(gv);
            cached_collision_positions_.push_back(pos);
            timesteps_to_collision_.push_back(std::numeric_limits<double>::infinity());
            particle_temperatures_.push_back(Tres);
            particle_omega_.push_back(phonon.mode_angular_frequency(mode));
            particle_occupation_.push_back(phonon.bose_occupation(Tres, mode));
            particle_energies_.push_back(0.0);
            particle_grid_id_.push_back(nearest_grid_index(geometry, pos));
            cached_collision_facets_.push_back(-1);
            cached_collision_conditions_.push_back('R');
            particle_alive_flags_.push_back(static_cast<std::uint8_t>(1));
            ++particle_count_;
            inserted.push_back({new_idx, dt_in[k]});
        }
    }
    return inserted;
}

// 函数说明：由粒子非平衡占据统计网格能量密度并做归一化修正。
void MonteCarloSolver::update_grid_energy_density(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    const int nsv = std::max(1, geometry.grid_count());
    grid_energy_density_.assign(static_cast<size_t>(nsv), 0.0);

#ifdef _OPENMP
    const int thread_count = std::max(1, omp_get_max_threads());
    ensure_tls_buffers(thread_count, nsv);
    std::fill(energy_tls_buffer_.begin(), energy_tls_buffer_.end(), 0.0);
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        double* local_energy = energy_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
#pragma omp for
        for (int i = 0; i < particle_count_; ++i) {
            const int sv = std::clamp(particle_grid_id_[i], 0, nsv - 1);
            const double n_eq = phonon.bose_occupation(grid_temperatures_[static_cast<size_t>(sv)], particle_modes_[i]);
            const double dn = particle_occupation_[i] - n_eq;
            particle_omega_[i] = phonon.mode_angular_frequency(particle_modes_[i]);
            particle_energies_[i] = 6.582119569e-4 * particle_omega_[i] * dn;  // hbar[eV*ps] * mode_angular_frequency[rad/ps] => eV
            local_energy[sv] += particle_energies_[i];
        }
    }
    for (int tid = 0; tid < thread_count; ++tid) {
        const double* local = energy_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
        for (int sv = 0; sv < nsv; ++sv) {
            grid_energy_density_[static_cast<size_t>(sv)] += local[sv];
        }
    }
#else
    for (int i = 0; i < particle_count_; ++i) {
        const int sv = std::clamp(particle_grid_id_[i], 0, nsv - 1);
        const double n_eq = phonon.bose_occupation(grid_temperatures_[static_cast<size_t>(sv)], particle_modes_[i]);
        const double dn = particle_occupation_[i] - n_eq;
        particle_omega_[i] = phonon.mode_angular_frequency(particle_modes_[i]);
        particle_energies_[i] = 6.582119569e-4 * particle_omega_[i] * dn;  // hbar[eV*ps] * mode_angular_frequency[rad/ps] => eV
        grid_energy_density_[static_cast<size_t>(sv)] += particle_energies_[i];
    }
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const int np = grid_particle_counts_[static_cast<size_t>(sv)];
        double norm_fac = 0.0;
        if (np > 0) {
            norm_fac = static_cast<double>(phonon.active_mode_count()) / static_cast<double>(np);
        }
        double e = grid_energy_density_[static_cast<size_t>(sv)] * norm_fac;
        e = phonon.normalize_to_energy_density(e);
        e += phonon.crystal_energy_density(grid_temperatures_[static_cast<size_t>(sv)]);
        grid_energy_density_[static_cast<size_t>(sv)] = e;
    }
}

// 函数说明：完成粒子计数、网格温度反演与粒子温度回写闭环。
void MonteCarloSolver::update_particle_temperatures(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    const int nsv = std::max(1, geometry.grid_count());
    grid_particle_counts_.assign(static_cast<size_t>(nsv), 0);
#ifdef _OPENMP
    const int thread_count = std::max(1, omp_get_max_threads());
    ensure_tls_buffers(thread_count, nsv);
    std::fill(count_tls_buffer_.begin(), count_tls_buffer_.end(), 0);
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        int* local_counts = count_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
#pragma omp for
        for (int i = 0; i < particle_count_; ++i) {
            const int sv = std::clamp(particle_grid_id_[i], 0, nsv - 1);
            local_counts[sv] += 1;
        }
    }
    for (int tid = 0; tid < thread_count; ++tid) {
        const int* local = count_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
        for (int sv = 0; sv < nsv; ++sv) {
            grid_particle_counts_[static_cast<size_t>(sv)] += local[sv];
        }
    }
#else
    for (int i = 0; i < particle_count_; ++i) {
        const int sv = std::clamp(particle_grid_id_[i], 0, nsv - 1);
        grid_particle_counts_[static_cast<size_t>(sv)] += 1;
    }
#endif

    update_grid_energy_density(geometry, phonon);
    apply_local_heat_source();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        grid_temperatures_[static_cast<size_t>(sv)] =
            phonon.temperature_from_energy_density(grid_energy_density_[static_cast<size_t>(sv)]);
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        const int sv = std::clamp(particle_grid_id_[i], 0, nsv - 1);
        particle_temperatures_[i] = grid_temperatures_[static_cast<size_t>(sv)];
    }
}

// 函数说明：按模态寿命推进粒子占据数，模拟弛豫散射过程。
void MonteCarloSolver::apply_lifetime_scattering(const PhononMaterial& phonon) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        const double n0 = phonon.bose_occupation(particle_temperatures_[i], particle_modes_[i]);
        const double tau = phonon.mode_lifetime(particle_temperatures_[i], particle_modes_[i]);
        if (tau > 0.0) {
            const double fac = std::exp(-time_step_ / tau);
            particle_occupation_[i] = n0 + (particle_occupation_[i] - n0) * fac;
        } else {
            particle_occupation_[i] = n0;
        }
    }
}

// 函数说明：统计网格热流并计算导热率估计（端点法与线性拟合法）。
void MonteCarloSolver::update_heat_flux_and_conductivity(const SimulationDomain& geometry) {
    const int nsv = std::max(1, geometry.grid_count());
    grid_heat_flux_.assign(static_cast<size_t>(nsv), {0.0, 0.0, 0.0});

#ifdef _OPENMP
    const int thread_count = std::max(1, omp_get_max_threads());
    ensure_tls_buffers(thread_count, nsv);
    std::fill(flux_tls_buffer_.begin(), flux_tls_buffer_.end(), Vec3 {0.0, 0.0, 0.0});
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        Vec3* local_flux = flux_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
#pragma omp for
        for (int i = 0; i < particle_count_; ++i) {
            const int sv = std::clamp(particle_grid_id_[i], 0, nsv - 1);
            local_flux[sv] = add(local_flux[sv], mul(particle_velocities_[i], particle_energies_[i]));
        }
    }
    for (int tid = 0; tid < thread_count; ++tid) {
        const Vec3* local = flux_tls_buffer_.data() + static_cast<size_t>(tid) * static_cast<size_t>(nsv);
        for (int sv = 0; sv < nsv; ++sv) {
            grid_heat_flux_[static_cast<size_t>(sv)] =
                add(grid_heat_flux_[static_cast<size_t>(sv)], local[sv]);
        }
    }
#else
    for (int i = 0; i < particle_count_; ++i) {
        const int sv = std::clamp(particle_grid_id_[i], 0, nsv - 1);
        grid_heat_flux_[static_cast<size_t>(sv)] =
            add(grid_heat_flux_[static_cast<size_t>(sv)], mul(particle_velocities_[i], particle_energies_[i]));
    }
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const int np = grid_particle_counts_[static_cast<size_t>(sv)];
        double norm_fac = 0.0;
        if (np > 0 && phonon_ != nullptr) {
            norm_fac = static_cast<double>(phonon_->active_mode_count()) / static_cast<double>(np);
        }
        grid_heat_flux_[static_cast<size_t>(sv)] = mul(grid_heat_flux_[static_cast<size_t>(sv)], norm_fac);
        if (phonon_ != nullptr) {
            grid_heat_flux_[static_cast<size_t>(sv)] = phonon_->normalize_to_energy_density(grid_heat_flux_[static_cast<size_t>(sv)]);
        }
        grid_heat_flux_[static_cast<size_t>(sv)] = mul(grid_heat_flux_[static_cast<size_t>(sv)], evpsa2_to_wm2_);
    }

    int axis = 0;
    if (!args_.grid_layout.empty() && args_.grid_layout.front() == "grid" && args_.grid_layout.size() >= 4) {
        const int nx = std::max(1, std::stoi(args_.grid_layout[1]));
        const int ny = std::max(1, std::stoi(args_.grid_layout[2]));
        const int nz = std::max(1, std::stoi(args_.grid_layout[3]));
        if (ny > nx && ny >= nz) {
            axis = 1;
        } else if (nz > nx && nz > ny) {
            axis = 2;
        }
    }

    const int total_np = std::max(1, std::accumulate(grid_particle_counts_.begin(), grid_particle_counts_.end(), 0));
    double phi_weighted = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:phi_weighted)
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        phi_weighted += grid_heat_flux_[static_cast<size_t>(sv)][axis] * static_cast<double>(grid_particle_counts_[static_cast<size_t>(sv)]);
    }
    average_heat_flux_along_axis_ = phi_weighted / static_cast<double>(total_np);

    if (!args_.compute_kappa) {
        thermal_conductivity_fit_ = 0.0;
        thermal_conductivity_endpoints_ = 0.0;
        thermal_conductivity_ = 0.0;
        return;
    }

    std::vector<double> T(static_cast<size_t>(nsv + 2), 0.0);
    for (int sv = 0; sv < nsv; ++sv) {
        T[static_cast<size_t>(sv + 1)] = grid_temperatures_[static_cast<size_t>(sv)];
    }
    double T_left = grid_temperatures_.front();
    double T_right = grid_temperatures_.back();
    const auto& rf = geometry.reservoir_facets();
    if (rf.size() >= 2) {
        const auto& fcent = geometry.mesh().facet_centroids();
        int f_left = rf.front();
        int f_right = rf.front();
        for (int f : rf) {
            if (f < 0 || f >= static_cast<int>(fcent.size())) {
                continue;
            }
            if (fcent[static_cast<size_t>(f)][axis] < fcent[static_cast<size_t>(f_left)][axis]) {
                f_left = f;
            }
            if (fcent[static_cast<size_t>(f)][axis] > fcent[static_cast<size_t>(f_right)][axis]) {
                f_right = f;
            }
        }
        T_left = geometry.reservoir_value_for_facet(f_left, T_left);
        T_right = geometry.reservoir_value_for_facet(f_right, T_right);
    }
    T.front() = T_left;
    T.back() = T_right;

    const auto& bmin = geometry.bounds_min();
    const auto& bmax = geometry.bounds_max();
    const double L = std::abs(bmax[axis] - bmin[axis]) * angstrom_to_meter_;

    // Method 1: endpoint temperature difference gradient.
    const double grad_end = (std::abs(L) > 1e-24) ? ((T_right - T_left) / L) : 0.0;
    if (std::abs(grad_end) > 1e-18) {
        thermal_conductivity_endpoints_ = -average_heat_flux_along_axis_ / grad_end;
    } else {
        thermal_conductivity_endpoints_ = 0.0;
    }

    // Method 2: linear fit on middle grids.
    const auto& centers = geometry.grid_centers();
    int i0 = 0;
    int i1 = nsv - 1;

    double sx = 0.0;
    double st = 0.0;
    double sxx = 0.0;
    double sxt = 0.0;
    int nfit = 0;
    for (int i = i0; i <= i1 && i < static_cast<int>(centers.size()); ++i) {
        const double x = centers[static_cast<size_t>(i)][axis] * angstrom_to_meter_;
        const double tt = grid_temperatures_[static_cast<size_t>(i)];
        sx += x;
        st += tt;
        sxx += x * x;
        sxt += x * tt;
        ++nfit;
    }
    double grad_fit = 0.0;
    const double dT_window = std::abs(
        grid_temperatures_[static_cast<size_t>(std::clamp(i1, 0, nsv - 1))] -
        grid_temperatures_[static_cast<size_t>(std::clamp(i0, 0, nsv - 1))]);
    if (nfit >= 2) {
        const double den = static_cast<double>(nfit) * sxx - sx * sx;
        if (std::abs(den) > 1e-30 && dT_window > 1e-6) {
            grad_fit = (static_cast<double>(nfit) * sxt - sx * st) / den;
        }
    }
    if (std::abs(grad_fit) > 1e-18) {
        thermal_conductivity_fit_ = -average_heat_flux_along_axis_ / grad_fit;
    } else {
        thermal_conductivity_fit_ = 0.0;
    }

    // Backward-compatible scalar uses endpoint gradient method.
    thermal_conductivity_ = thermal_conductivity_endpoints_;
}

// 函数说明：在单个时间步内推进粒子运动并处理可能的多次边界碰撞。
void MonteCarloSolver::advance_particle(const SimulationDomain& geometry, const PhononMaterial& phonon, int i, double dt_remaining) {
    if (i < 0 || i >= particle_count_) {
        return;
    }
    if (i < static_cast<int>(particle_alive_flags_.size()) && particle_alive_flags_[static_cast<size_t>(i)] == 0) {
        return;
    }
    double remaining = dt_remaining;
    int guard = 0;
    while (remaining > 1e-14 && guard < 64) {
        ++guard;
        const double t_hit = std::isinf(timesteps_to_collision_[i]) ? std::numeric_limits<double>::infinity() : timesteps_to_collision_[i] * time_step_;
        if (!std::isfinite(t_hit) || t_hit > remaining) {
            particle_positions_[i] = add(particle_positions_[i], mul(particle_velocities_[i], remaining));
            if (std::isfinite(timesteps_to_collision_[i])) {
                timesteps_to_collision_[i] -= remaining / time_step_;
            }
            remaining = 0.0;
        } else {
            // Move to collision and process boundary event.
            particle_positions_[i] = cached_collision_positions_[i];
            remaining -= std::max(0.0, t_hit);
            process_boundary_collision(geometry, i);
            if (particle_alive_flags_[static_cast<size_t>(i)] == 0) {
                break;
            }
            update_collision_cache_single(geometry, i);
            if (timesteps_to_collision_[i] * time_step_ < 1e-12) {
                particle_positions_[i] = add(particle_positions_[i], mul(particle_velocities_[i], 1e-12));
                update_collision_cache_single(geometry, i);
            }
        }
    }
    (void) phonon;
}

// 函数说明：初始化收敛输出文件表头，定义温度与热输运列。
void MonteCarloSolver::write_convergence_header() {
    if (args_.output_folder.empty()) {
        return;
    }
    std::ofstream out(std::filesystem::path(args_.output_folder) / "convergence.txt", std::ios::trunc);
    out << "# timestep time_ps"; 
    const int ngrid = static_cast<int>(grid_temperatures_.size());
    for (int i = 0; i < ngrid; ++i) {
        out << " T_" << (i + 1);
    }
    out << " heatflux kappa_fit kappa_end\n";
}

// 函数说明：追加当前时间步温度、热流与导热率统计结果。
void MonteCarloSolver::append_convergence_row() const {
    if (args_.output_folder.empty()) {
        return;
    }
    std::ofstream out(std::filesystem::path(args_.output_folder) / "convergence.txt", std::ios::app);
    // 修改点：在第一列后增加 elapsed_time_
    out << current_timestep_ << " " << elapsed_time_; 
    for (double tsv : grid_temperatures_) {
        out << " " << tsv;
    }
    out << " " << average_heat_flux_along_axis_ << " " << thermal_conductivity_fit_ << " " << thermal_conductivity_endpoints_ << '\n';
}

void MonteCarloSolver::append_profile_summary(std::ostream& out) const {
    out << "\n[performance]\n";
#ifdef _OPENMP
    out << "openmp_enabled = true\n";
#else
    out << "openmp_enabled = false\n";
#endif
    out << "openmp_thread_count = " << openmp_thread_count_ << '\n';
    out << "profile_timers_enabled = " << (profile_timers_enabled_ ? "true" : "false") << '\n';
    const long long rough_total = rough_events_total_.load(std::memory_order_relaxed);
    const long long rough_spec = rough_specular_selected_.load(std::memory_order_relaxed);
    const long long rough_diff = rough_diffuse_selected_.load(std::memory_order_relaxed);
    const long long rough_fb_no_data = rough_fallback_missing_rough_data_.load(std::memory_order_relaxed);
    const long long rough_fb_no_spec = rough_fallback_missing_spec_match_.load(std::memory_order_relaxed);
    const long long rough_fb_pool = rough_fallback_outgoing_pool_.load(std::memory_order_relaxed);
    const long long rough_fb_rand = rough_fallback_global_random_.load(std::memory_order_relaxed);
    const double rough_den = std::max(1.0, static_cast<double>(rough_total));
    out << "\n[rough_boundary_diagnostics]\n";
    out << "events_total = " << rough_total << '\n';
    out << "specular_selected = " << rough_spec << '\n';
    out << "diffuse_selected = " << rough_diff << '\n';
    out << "specular_ratio = " << (static_cast<double>(rough_spec) / rough_den) << '\n';
    out << "diffuse_ratio = " << (static_cast<double>(rough_diff) / rough_den) << '\n';
    out << "fallback_missing_rough_data = " << rough_fb_no_data << '\n';
    out << "fallback_missing_spec_match = " << rough_fb_no_spec << '\n';
    out << "fallback_outgoing_pool = " << rough_fb_pool << '\n';
    out << "fallback_global_random = " << rough_fb_rand << '\n';
    out << "fallback_total = " << (rough_fb_no_data + rough_fb_no_spec + rough_fb_pool + rough_fb_rand) << '\n';
    out << "mode_map_csv = disabled\n";
    out << "\n[escaped_particle_recovery]\n";
    out << "check_interval = " << escaped_recovery_check_interval_ << '\n';
    out << "recovery_events = " << escaped_recovery_events_ << '\n';
    out << "recovered_particles_total = " << escaped_recovered_particles_ << '\n';
    if (!profile_timers_enabled_) {
        return;
    }
    out << "profile_total_seconds = " << timer_total_ << '\n';
    const double denom = std::max(timer_total_, 1e-18);
    auto write_seg = [&out, denom](const char* name, double sec) {
        out << name << "_seconds = " << sec << '\n';
        out << name << "_percent = " << (100.0 * sec / denom) << '\n';
    };
    write_seg("advance_main", timer_advance_main_);
    write_seg("remove_absorb_1", timer_remove_absorb_1_);
    write_seg("inject_build", timer_inject_build_);
    write_seg("inject_cache", timer_inject_cache_);
    write_seg("advance_injected", timer_advance_injected_);
    write_seg("remove_absorb_2", timer_remove_absorb_2_);
    write_seg("update_temp", timer_update_temp_);
    write_seg("lifetime", timer_lifetime_);
    write_seg("stats", timer_stats_);
}

void MonteCarloSolver::report_timestep_timers_if_needed() const {
    if (!profile_timers_enabled_ || current_timestep_ <= 0) {
        return;
    }
    if (current_timestep_ % 100 != 0 && current_timestep_ != args_.iterations) {
        return;
    }
    const double denom = std::max(timer_total_, 1e-18);
    auto pct = [denom](double x) { return 100.0 * x / denom; };
    std::cout << "[profile] timesteps=" << current_timestep_
              << " total=" << timer_total_ << "s\n";
    std::cout << "  advance_main=" << pct(timer_advance_main_) << "%\n";
    std::cout << "  remove_absorb_1=" << pct(timer_remove_absorb_1_) << "%\n";
    std::cout << "  inject_build=" << pct(timer_inject_build_) << "%\n";
    std::cout << "  inject_cache=" << pct(timer_inject_cache_) << "%\n";
    std::cout << "  advance_injected=" << pct(timer_advance_injected_) << "%\n";
    std::cout << "  remove_absorb_2=" << pct(timer_remove_absorb_2_) << "%\n";
    std::cout << "  update_temp=" << pct(timer_update_temp_) << "%\n";
    std::cout << "  lifetime=" << pct(timer_lifetime_) << "%\n";
    std::cout << "  stats=" << pct(timer_stats_) << "%\n";
}

// 函数说明：执行一次完整时间步流程：推进、注入、温度更新、散射与输出。
void MonteCarloSolver::run_timestep() {
    if (geometry_ == nullptr) {
        throw std::runtime_error("MonteCarloSolver geometry is not set.");
    }
    if (phonon_ == nullptr) {
        throw std::runtime_error("MonteCarloSolver phonon is not set.");
    }
    const SimulationDomain& geometry = *geometry_;
    const PhononMaterial& phonon = *phonon_;
    using Clock = std::chrono::steady_clock;
    const auto t_step_begin = Clock::now();

    const int n_before = particle_count_;
    const auto t_advance_begin = Clock::now();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int i = 0; i < n_before; ++i) {
        advance_particle(geometry, phonon, i, time_step_);
        if (i < particle_count_ && i < static_cast<int>(particle_alive_flags_.size()) && particle_alive_flags_[static_cast<size_t>(i)] != 0) {
            particle_grid_id_[i] = nearest_grid_index(geometry, particle_positions_[i]);
        }
    }
    const auto t_advance_end = Clock::now();
    const auto t_remove1_begin = t_advance_end;
    remove_absorbed_particles();
    const auto t_remove1_end = Clock::now();

    const auto t_inject_begin = t_remove1_end;
    const auto injected = inject_particles_from_reservoirs(geometry, phonon);
    const auto t_inject_end = Clock::now();
    if (!injected.empty()) {
        const auto t_cache_begin = Clock::now();
        std::vector<int> new_idx;
        new_idx.reserve(injected.size());
        for (const auto& [idx, _] : injected) {
            new_idx.push_back(idx);
        }
        update_collision_cache(geometry, new_idx);
        const auto t_cache_end = Clock::now();
        const auto t_adv_inj_begin = t_cache_end;
        const int ninj = static_cast<int>(injected.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (int k = 0; k < ninj; ++k) {
            const auto [idx, dt_in] = injected[static_cast<size_t>(k)];
            const double remain = std::max(0.0, time_step_ - dt_in);
            if (remain > 1e-14 && idx < particle_count_) {
                advance_particle(geometry, phonon, idx, remain);
            }
            if (idx < particle_count_ &&
                idx < static_cast<int>(particle_alive_flags_.size()) &&
                particle_alive_flags_[static_cast<size_t>(idx)] != 0) {
                particle_grid_id_[idx] = nearest_grid_index(geometry, particle_positions_[idx]);
            }
        }
        const auto t_adv_inj_end = Clock::now();
        const auto t_remove2_begin = t_adv_inj_end;
        remove_absorbed_particles();
        const auto t_remove2_end = Clock::now();
        if (profile_timers_enabled_) {
            timer_inject_cache_ += std::chrono::duration<double>(t_cache_end - t_cache_begin).count();
            timer_advance_injected_ += std::chrono::duration<double>(t_adv_inj_end - t_adv_inj_begin).count();
            timer_remove_absorb_2_ += std::chrono::duration<double>(t_remove2_end - t_remove2_begin).count();
        }
    }

    if (escaped_recovery_check_interval_ > 0 &&
        (current_timestep_ % escaped_recovery_check_interval_) == 0) {
        recover_escaped_particles(geometry);
    }

    const auto t_temp_begin = Clock::now();
    update_particle_temperatures(geometry, phonon);
    const auto t_temp_end = Clock::now();
    const auto t_life_begin = t_temp_end;
    apply_lifetime_scattering(phonon);
    const auto t_life_end = Clock::now();

    ++current_timestep_;
    elapsed_time_ += time_step_;
    const auto t_stats_begin = Clock::now();
    if (current_timestep_ % convergence_write_interval_ == 0 || current_timestep_ == 1) {
        update_heat_flux_and_conductivity(geometry);
        append_convergence_row();
    }
    const auto t_stats_end = Clock::now();
    const auto t_step_end = t_stats_end;
    if (profile_timers_enabled_) {
        timer_total_ += std::chrono::duration<double>(t_step_end - t_step_begin).count();
        timer_advance_main_ += std::chrono::duration<double>(t_advance_end - t_advance_begin).count();
        timer_remove_absorb_1_ += std::chrono::duration<double>(t_remove1_end - t_remove1_begin).count();
        timer_inject_build_ += std::chrono::duration<double>(t_inject_end - t_inject_begin).count();
        timer_update_temp_ += std::chrono::duration<double>(t_temp_end - t_temp_begin).count();
        timer_lifetime_ += std::chrono::duration<double>(t_life_end - t_life_begin).count();
        timer_stats_ += std::chrono::duration<double>(t_stats_end - t_stats_begin).count();
    }
    int total_iters = args_.iterations;
    int print_interval = std::max(1, total_iters / 100);

    if (current_timestep_ % print_interval == 0 || current_timestep_ == total_iters) {
        double progress = (static_cast<double>(current_timestep_) / total_iters) * 100.0;
        
        std::cout << "--- Progress: " << std::fixed << std::setprecision(1) << progress << "% ---" << std::endl;
        
        std::cout << "Temperature Profile (K): ";
        for (double t : grid_temperatures_) {
            std::cout << std::setprecision(2) << t << " ";
        }
        std::cout << std::endl;

        if (args_.compute_kappa) {
            std::cout << "Current Conductivity (Fit): " << thermal_conductivity_fit_ << " W/mK" << std::endl;
            std::cout << "Current Conductivity (End): " << thermal_conductivity_endpoints_ << " W/mK" << std::endl;
        }
        std::cout << "-------------------------" << std::endl;
    }
    report_timestep_timers_if_needed();
}
