#include "MonteCarloSolver.h"

#include "SimulationDomain.h"
#include "PhononMaterial.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

MonteCarloSolver::Vec3 MonteCarloSolver::add(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
MonteCarloSolver::Vec3 MonteCarloSolver::sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
MonteCarloSolver::Vec3 MonteCarloSolver::mul(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double MonteCarloSolver::dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double MonteCarloSolver::norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

MonteCarloSolver::MonteCarloSolver(const SimulationConfig& args, const SimulationDomain& geometry, const PhononMaterial& phonon)
    : args_(args), geometry_(&geometry), phonon_(&phonon) {
    particle_count_ = std::max(1, static_cast<int>(std::llround(args_.particle_count)));
    time_step_ = std::max(1e-12, args_.time_step);
    push_eps_ = 1e-10 * std::max(time_step_, 1.0);
    particle_density_ = static_cast<double>(particle_count_) / std::max(geometry.volume(), 1e-12);

    std::cout << "MonteCarloSolver initialized: particle_count=" << particle_count_
              << ", time_step=" << time_step_
              << ", density=" << particle_density_ << '\n';
    std::cout << "Thermal conductivity estimation: "
              << (args_.compute_thermal_conductivity ? "enabled" : "disabled")
              << '\n';
#ifdef _OPENMP
    std::cout << "OpenMP enabled: max_threads=" << omp_get_max_threads() << '\n';
#else
    std::cout << "OpenMP enabled: no\n";
#endif

    if (!args_.results_base_folder.empty()) {
        std::filesystem::create_directories(args_.results_base_folder);
    }

    initialize_particles(geometry, phonon);
    initialize_local_heat_source(geometry);
    write_convergence_header();
    update_heat_flux_and_conductivity(geometry);
    append_convergence_row();
}

void MonteCarloSolver::initialize_particles(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    const auto& mesh = geometry.mesh();
    particle_positions_ = mesh.sample_volume_points(particle_count_);
    particle_subvolume_id_.assign(static_cast<size_t>(particle_count_), 0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        particle_subvolume_id_[i] = nearest_subvolume_index(geometry, particle_positions_[i]);
    }

    initialize_particle_modes(phonon);
    initialize_particle_temperatures(geometry);
    initialize_particle_velocities(phonon);
    initialize_reservoir_injection(geometry, phonon);
    initialize_rough_boundary_scattering(geometry, phonon);
    particle_omega_.resize(static_cast<size_t>(particle_count_));
    particle_occupation_.resize(static_cast<size_t>(particle_count_));
    particle_energies_.resize(static_cast<size_t>(particle_count_));
    particle_alive_flags_.assign(static_cast<size_t>(particle_count_), static_cast<std::uint8_t>(1));
    const double tmin = particle_temperatures_.empty() ? 300.0 : *std::min_element(particle_temperatures_.begin(), particle_temperatures_.end());
    const int nsv = std::max(1, geometry.subvolume_count());
    subvolume_temperatures_.assign(static_cast<size_t>(nsv), tmin);
    subvolume_particle_counts_.assign(static_cast<size_t>(nsv), 0);
    subvolume_energy_density_.assign(static_cast<size_t>(nsv), 0.0);
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
    cached_collision_conditions_.assign(static_cast<size_t>(particle_count_), 'P');
    std::vector<int> all_idx(static_cast<size_t>(particle_count_));
    std::iota(all_idx.begin(), all_idx.end(), 0);
    update_collision_cache(geometry, all_idx);
    update_particle_temperatures(geometry, phonon);
}

void MonteCarloSolver::initialize_particle_modes(const PhononMaterial& phonon) {
    particle_modes_.resize(static_cast<size_t>(particle_count_));
    for (int i = 0; i < particle_count_; ++i) {
        particle_modes_[i] = phonon.sample_active_mode(rng_);
    }
}

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
    if (!args_.initial_temperature_profile.empty()) {
        key = args_.initial_temperature_profile.front();
    }

    if (key == "cold") {
        std::fill(particle_temperatures_.begin(), particle_temperatures_.end(), tmin);
    } else if (key == "hot") {
        std::fill(particle_temperatures_.begin(), particle_temperatures_.end(), tmax);
    } else if (key == "mean") {
        std::fill(particle_temperatures_.begin(), particle_temperatures_.end(), tmean);
    } else if (key == "random") {
        std::uniform_real_distribution<double> U(tmin, tmax);
        for (auto& t : particle_temperatures_) {
            t = U(rng_);
        }
    } else if (key == "linear") {
        const auto& centers = geometry.subvolume_centers();
        int axis = 0;
        if (!args_.subvolume_layout.empty() && args_.subvolume_layout.front() == "slice" && args_.subvolume_layout.size() >= 3) {
            axis = std::clamp(std::stoi(args_.subvolume_layout[2]), 0, 2);
        } else if (!centers.empty()) {
            Vec3 mn = centers.front();
            Vec3 mx = centers.front();
            for (const auto& c : centers) {
                for (int k = 0; k < 3; ++k) {
                    mn[k] = std::min(mn[k], c[k]);
                    mx[k] = std::max(mx[k], c[k]);
                }
            }
            const Vec3 ex = sub(mx, mn);
            axis = static_cast<int>(std::max_element(ex.begin(), ex.end()) - ex.begin());
        }
        double pmin = particle_positions_.front()[axis];
        double pmax = particle_positions_.front()[axis];
        for (const auto& p : particle_positions_) {
            pmin = std::min(pmin, p[axis]);
            pmax = std::max(pmax, p[axis]);
        }
        const double den = std::max(1e-12, pmax - pmin);
        for (int i = 0; i < particle_count_; ++i) {
            const double s = (particle_positions_[i][axis] - pmin) / den;
            particle_temperatures_[i] = tmin + (tmax - tmin) * s;
        }
    } else {
        std::fill(particle_temperatures_.begin(), particle_temperatures_.end(), tmean);
    }
}

MonteCarloSolver::Vec3 MonteCarloSolver::random_unit_vector() {
    std::normal_distribution<double> N(0.0, 1.0);
    Vec3 v {N(rng_), N(rng_), N(rng_)};
    const double n = norm(v);
    if (n <= 1e-12) {
        return {1.0, 0.0, 0.0};
    }
    return mul(v, 1.0 / n);
}

void MonteCarloSolver::initialize_particle_velocities(const PhononMaterial& phonon) {
    particle_velocities_.resize(static_cast<size_t>(particle_count_));
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        particle_velocities_[i] = phonon.mode_group_velocity(particle_modes_[i]);
    }
}

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
            reservoir_emit_counter_[static_cast<size_t>(r)][m] = U(rng_);
        }
    }
}

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

    for (int facet = 0; facet < nfacets; ++facet) {
        if (!geometry.is_rough_facet(facet)) {
            continue;
        }
        RoughFacetData rd;
        rd.facet = facet;
        rd.specularity.assign(static_cast<size_t>(na), 0.0);
        rd.spec_match_active.assign(static_cast<size_t>(na), -1);
        rd.diffuse_begin.assign(static_cast<size_t>(na), -1);
        rd.diffuse_end.assign(static_cast<size_t>(na), -1);

        const double eta = std::max(0.0, geometry.roughness_for_facet(facet, 0.0));
        const Vec3 n_out = mesh.facet_normals()[static_cast<size_t>(facet)];
        const Vec3 n_in {-n_out[0], -n_out[1], -n_out[2]};
        std::vector<int> incoming;
        std::vector<int> outgoing;
        incoming.reserve(static_cast<size_t>(na));
        outgoing.reserve(static_cast<size_t>(na));

        for (int ai = 0; ai < na; ++ai) {
            const Vec3 v = v_active[static_cast<size_t>(ai)];
            const double vn = dot(v, n_in);
            if (vn < 0.0) {
                incoming.push_back(ai);
            } else if (vn > 0.0) {
                outgoing.push_back(ai);
            }
            const double vnorm = std::max(1e-12, vnorm_active[static_cast<size_t>(ai)]);
            const double incidence_cos = std::abs(vn) / vnorm;
            const double x = 2.0 * eta * incidence_cos * knorm_active[static_cast<size_t>(ai)];
            double p = std::exp(-(x * x));
            if (!std::isfinite(p)) {
                p = 0.0;
            }
            rd.specularity[static_cast<size_t>(ai)] = std::clamp(p, 0.0, 1.0);
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

        const int rid = static_cast<int>(rough_boundary_data_.size());
        rough_boundary_data_.push_back(std::move(rd));
        facet_to_rough_data_[static_cast<size_t>(facet)] = rid;
    }
}

void MonteCarloSolver::initialize_local_heat_source(const SimulationDomain& geometry) {
    local_heat_source_enabled_ = false;
    local_heat_source_subvolume_mask_.clear();
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

    const std::string mode = args_.heat_source_mode;
    if (mode == "relative") {
        const auto& bmin = geometry.bounds_min();
        const auto& bmax = geometry.bounds_max();
        const Vec3 ext {bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]};
        for (int k = 0; k < 3; ++k) {
            local_heat_source_min_[k] = bmin[k] + local_heat_source_min_[k] * ext[k];
            local_heat_source_max_[k] = bmin[k] + local_heat_source_max_[k] * ext[k];
        }
    } else if (mode != "absolute") {
        throw std::runtime_error("Unsupported heat_source_mode. Use relative or absolute.");
    }

    for (int k = 0; k < 3; ++k) {
        if (local_heat_source_min_[k] > local_heat_source_max_[k]) {
            std::swap(local_heat_source_min_[k], local_heat_source_max_[k]);
        }
    }

    const auto& centers = geometry.subvolume_centers();
    local_heat_source_subvolume_mask_.assign(centers.size(), static_cast<std::uint8_t>(0));
    int selected = 0;
    for (size_t i = 0; i < centers.size(); ++i) {
        const Vec3 c = centers[i];
        const bool inside =
            c[0] >= local_heat_source_min_[0] && c[0] <= local_heat_source_max_[0] &&
            c[1] >= local_heat_source_min_[1] && c[1] <= local_heat_source_max_[1] &&
            c[2] >= local_heat_source_min_[2] && c[2] <= local_heat_source_max_[2];
        if (inside) {
            local_heat_source_subvolume_mask_[i] = static_cast<std::uint8_t>(1);
            ++selected;
        }
    }

    if (selected == 0) {
        std::cerr << "Warning: heat source region does not include any subvolume center. Ignoring local heat source.\n";
        return;
    }

    local_heat_source_enabled_ = true;
    std::cout << "Local heat source enabled: selected_subvolumes=" << selected
              << ", power_density=" << local_heat_source_power_density_wm3_ << " W/m^3\n";
}

void MonteCarloSolver::apply_local_heat_source() {
    if (!local_heat_source_enabled_) {
        return;
    }
    const double delta_e = local_heat_source_power_density_wm3_ * wm3_to_evpsa3_ * time_step_;
    if (!std::isfinite(delta_e) || std::abs(delta_e) <= 0.0) {
        return;
    }
    const size_t n = std::min(subvolume_energy_density_.size(), local_heat_source_subvolume_mask_.size());
    for (size_t i = 0; i < n; ++i) {
        if (local_heat_source_subvolume_mask_[i] != 0) {
            subvolume_energy_density_[i] += delta_e;
        }
    }
}

int MonteCarloSolver::sample_diffuse_active_mode(int rough_idx, int in_ai) const {
    const int na = (phonon_ != nullptr) ? phonon_->active_mode_count() : 0;
    if (na <= 0) {
        return 0;
    }
    if (rough_idx < 0 || rough_idx >= static_cast<int>(rough_boundary_data_.size())) {
        std::uniform_int_distribution<int> U(0, na - 1);
        return U(rng_);
    }
    const auto& rd = rough_boundary_data_[static_cast<size_t>(rough_idx)];
    if (in_ai >= 0 &&
        in_ai < static_cast<int>(rd.diffuse_begin.size()) &&
        in_ai < static_cast<int>(rd.diffuse_end.size()) &&
        !rd.outgoing_sorted_active.empty()) {
        const int ib = rd.diffuse_begin[static_cast<size_t>(in_ai)];
        const int ie = rd.diffuse_end[static_cast<size_t>(in_ai)];
        if (ie > ib && ib >= 0 && ie <= static_cast<int>(rd.outgoing_sorted_active.size())) {
            std::uniform_int_distribution<int> U(ib, ie - 1);
            return rd.outgoing_sorted_active[static_cast<size_t>(U(rng_))];
        }
    }
    if (!rd.outgoing_active.empty()) {
        std::uniform_int_distribution<int> U(0, static_cast<int>(rd.outgoing_active.size()) - 1);
        return rd.outgoing_active[static_cast<size_t>(U(rng_))];
    }
    std::uniform_int_distribution<int> U(0, na - 1);
    return U(rng_);
}

std::array<int, 2> MonteCarloSolver::select_reflected_mode(
    const SimulationDomain& geometry,
    const PhononMaterial& phonon,
    int rough_idx,
    const std::array<int, 2>& in_mode,
    const Vec3& collision_pos,
    double& out_occupation,
    double in_occupation) const {
    const int na = phonon.active_mode_count();
    if (rough_idx < 0 || rough_idx >= static_cast<int>(rough_boundary_data_.size()) || na <= 0) {
        out_occupation = in_occupation;
        return in_mode;
    }
    const auto& rd = rough_boundary_data_[static_cast<size_t>(rough_idx)];
    int in_ai = phonon.active_index_for_mode(in_mode);
    if (in_ai < 0 || in_ai >= na) {
        in_ai = -1;
    }

    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double p_spec = (in_ai >= 0 && in_ai < static_cast<int>(rd.specularity.size()))
        ? std::clamp(rd.specularity[static_cast<size_t>(in_ai)], 0.0, 1.0)
        : 0.0;
    if (in_ai >= 0 && U01(rng_) <= p_spec) {
        const int out_ai = rd.spec_match_active[static_cast<size_t>(in_ai)];
        if (out_ai >= 0 && out_ai < na) {
            out_occupation = in_occupation;
            return phonon.active_mode_at(out_ai);
        }
    }

    const int out_ai = sample_diffuse_active_mode(rough_idx, in_ai);
    const std::array<int, 2> out_mode = phonon.active_mode_at(out_ai);
    const int nsv = std::max(1, geometry.subvolume_count());
    const int sv = std::clamp(nearest_subvolume_index(geometry, collision_pos), 0, nsv - 1);
    const double Tdiff = (sv >= 0 && sv < static_cast<int>(subvolume_temperatures_.size()))
        ? subvolume_temperatures_[static_cast<size_t>(sv)] : 300.0;
    out_occupation = phonon.bose_occupation(Tdiff, out_mode);
    return out_mode;
}

int MonteCarloSolver::nearest_subvolume_index(const SimulationDomain& geometry, const Vec3& p) const {
    const auto& centers = geometry.subvolume_centers();
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

void MonteCarloSolver::update_collision_cache(const SimulationDomain& geometry, const std::vector<int>& indices) {
    const auto& mesh = geometry.mesh();
    const int nidx = static_cast<int>(indices.size());
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int k = 0; k < nidx; ++k) {
        const int i = indices[static_cast<size_t>(k)];
        if (i < 0 || i >= particle_count_) {
            continue;
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
}

void MonteCarloSolver::process_boundary_collision(const SimulationDomain& geometry, int i) {
    const int facet = cached_collision_facets_[i];
    const char cond = cached_collision_conditions_[i];

    particle_positions_[i] = cached_collision_positions_[i];
    if (facet < 0 || facet >= static_cast<int>(geometry.mesh().facet_normals().size())) {
        return;
    }

    const Vec3 n = geometry.mesh().facet_normals()[facet];
    const double vn = dot(particle_velocities_[i], n);

    if (cond == 'T' || cond == 'F') {
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
            const double p_spec = compute_roughness_specularity(geometry, *phonon_, i, facet);
            std::uniform_real_distribution<double> U(0.0, 1.0);
            if (U(rng_) <= p_spec) {
                particle_velocities_[i] = sub(particle_velocities_[i], mul(n, 2.0 * vn));
            } else {
                particle_modes_[i] = phonon_->sample_active_mode(rng_);
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
    particle_subvolume_id_[i] = nearest_subvolume_index(geometry, particle_positions_[i]);
}

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
    remap_veci(particle_subvolume_id_);
    remap_veci(cached_collision_facets_);
    remap_vecc(cached_collision_conditions_);
    particle_alive_flags_.assign(alive_count, static_cast<std::uint8_t>(1));
    particle_count_ = static_cast<int>(alive_count);
}

std::vector<std::pair<int, double>> MonteCarloSolver::inject_particles_from_reservoirs(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    std::vector<std::pair<int, double>> inserted;
    if (reservoir_count_ <= 0 || reservoir_modes_.empty()) {
        return inserted;
    }

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
                dt_in.push_back(time_step_ * U01(rng_));
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
            particle_subvolume_id_.push_back(nearest_subvolume_index(geometry, pos));
            cached_collision_facets_.push_back(-1);
            cached_collision_conditions_.push_back('P');
            particle_alive_flags_.push_back(static_cast<std::uint8_t>(1));
            ++particle_count_;
            inserted.push_back({new_idx, dt_in[k]});
        }
    }
    return inserted;
}

void MonteCarloSolver::update_subvolume_energy_density(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    const int nsv = std::max(1, geometry.subvolume_count());
    subvolume_energy_density_.assign(static_cast<size_t>(nsv), 0.0);

#ifdef _OPENMP
    const int thread_count = omp_get_max_threads();
    std::vector<std::vector<double>> energy_tls(
        static_cast<size_t>(thread_count),
        std::vector<double>(static_cast<size_t>(nsv), 0.0));
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& local_energy = energy_tls[static_cast<size_t>(tid)];
#pragma omp for
        for (int i = 0; i < particle_count_; ++i) {
            const int sv = std::clamp(particle_subvolume_id_[i], 0, nsv - 1);
            const double n_eq = phonon.bose_occupation(subvolume_temperatures_[static_cast<size_t>(sv)], particle_modes_[i]);
            const double dn = particle_occupation_[i] - n_eq;
            particle_omega_[i] = phonon.mode_angular_frequency(particle_modes_[i]);
            particle_energies_[i] = 6.582119569e-4 * particle_omega_[i] * dn;  // hbar[eV*ps] * mode_angular_frequency[rad/ps] => eV
            local_energy[static_cast<size_t>(sv)] += particle_energies_[i];
        }
    }
    for (const auto& local : energy_tls) {
        for (int sv = 0; sv < nsv; ++sv) {
            subvolume_energy_density_[static_cast<size_t>(sv)] += local[static_cast<size_t>(sv)];
        }
    }
#else
    for (int i = 0; i < particle_count_; ++i) {
        const int sv = std::clamp(particle_subvolume_id_[i], 0, nsv - 1);
        const double n_eq = phonon.bose_occupation(subvolume_temperatures_[static_cast<size_t>(sv)], particle_modes_[i]);
        const double dn = particle_occupation_[i] - n_eq;
        particle_omega_[i] = phonon.mode_angular_frequency(particle_modes_[i]);
        particle_energies_[i] = 6.582119569e-4 * particle_omega_[i] * dn;  // hbar[eV*ps] * mode_angular_frequency[rad/ps] => eV
        subvolume_energy_density_[static_cast<size_t>(sv)] += particle_energies_[i];
    }
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const int np = subvolume_particle_counts_[static_cast<size_t>(sv)];
        double norm_fac = 0.0;
        if (np > 0) {
            norm_fac = static_cast<double>(phonon.active_mode_count()) / static_cast<double>(np);
        }
        double e = subvolume_energy_density_[static_cast<size_t>(sv)] * norm_fac;
        e = phonon.normalize_to_energy_density(e);
        e += phonon.crystal_energy_density(subvolume_temperatures_[static_cast<size_t>(sv)]);
        subvolume_energy_density_[static_cast<size_t>(sv)] = e;
    }
}

void MonteCarloSolver::update_particle_temperatures(const SimulationDomain& geometry, const PhononMaterial& phonon) {
    const int nsv = std::max(1, geometry.subvolume_count());
    subvolume_particle_counts_.assign(static_cast<size_t>(nsv), 0);
#ifdef _OPENMP
    const int thread_count = omp_get_max_threads();
    std::vector<std::vector<int>> count_tls(
        static_cast<size_t>(thread_count),
        std::vector<int>(static_cast<size_t>(nsv), 0));
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& local_counts = count_tls[static_cast<size_t>(tid)];
#pragma omp for
        for (int i = 0; i < particle_count_; ++i) {
            particle_subvolume_id_[i] = nearest_subvolume_index(geometry, particle_positions_[i]);
            const int sv = std::clamp(particle_subvolume_id_[i], 0, nsv - 1);
            local_counts[static_cast<size_t>(sv)] += 1;
        }
    }
    for (const auto& local : count_tls) {
        for (int sv = 0; sv < nsv; ++sv) {
            subvolume_particle_counts_[static_cast<size_t>(sv)] += local[static_cast<size_t>(sv)];
        }
    }
#else
    for (int i = 0; i < particle_count_; ++i) {
        particle_subvolume_id_[i] = nearest_subvolume_index(geometry, particle_positions_[i]);
        const int sv = std::clamp(particle_subvolume_id_[i], 0, nsv - 1);
        subvolume_particle_counts_[static_cast<size_t>(sv)] += 1;
    }
#endif

    update_subvolume_energy_density(geometry, phonon);
    apply_local_heat_source();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        subvolume_temperatures_[static_cast<size_t>(sv)] =
            phonon.temperature_from_energy_density(subvolume_energy_density_[static_cast<size_t>(sv)]);
    }
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < particle_count_; ++i) {
        const int sv = std::clamp(particle_subvolume_id_[i], 0, nsv - 1);
        particle_temperatures_[i] = subvolume_temperatures_[static_cast<size_t>(sv)];
    }
}

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

void MonteCarloSolver::update_heat_flux_and_conductivity(const SimulationDomain& geometry) {
    const int nsv = std::max(1, geometry.subvolume_count());
    subvolume_heat_flux_.assign(static_cast<size_t>(nsv), {0.0, 0.0, 0.0});

#ifdef _OPENMP
    const int thread_count = omp_get_max_threads();
    std::vector<std::vector<Vec3>> flux_tls(
        static_cast<size_t>(thread_count),
        std::vector<Vec3>(static_cast<size_t>(nsv), Vec3 {0.0, 0.0, 0.0}));
#pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        auto& local_flux = flux_tls[static_cast<size_t>(tid)];
#pragma omp for
        for (int i = 0; i < particle_count_; ++i) {
            const int sv = std::clamp(particle_subvolume_id_[i], 0, nsv - 1);
            local_flux[static_cast<size_t>(sv)] =
                add(local_flux[static_cast<size_t>(sv)], mul(particle_velocities_[i], particle_energies_[i]));
        }
    }
    for (const auto& local : flux_tls) {
        for (int sv = 0; sv < nsv; ++sv) {
            subvolume_heat_flux_[static_cast<size_t>(sv)] =
                add(subvolume_heat_flux_[static_cast<size_t>(sv)], local[static_cast<size_t>(sv)]);
        }
    }
#else
    for (int i = 0; i < particle_count_; ++i) {
        const int sv = std::clamp(particle_subvolume_id_[i], 0, nsv - 1);
        subvolume_heat_flux_[static_cast<size_t>(sv)] =
            add(subvolume_heat_flux_[static_cast<size_t>(sv)], mul(particle_velocities_[i], particle_energies_[i]));
    }
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        const int np = subvolume_particle_counts_[static_cast<size_t>(sv)];
        double norm_fac = 0.0;
        if (np > 0 && phonon_ != nullptr) {
            norm_fac = static_cast<double>(phonon_->active_mode_count()) / static_cast<double>(np);
        }
        subvolume_heat_flux_[static_cast<size_t>(sv)] = mul(subvolume_heat_flux_[static_cast<size_t>(sv)], norm_fac);
        if (phonon_ != nullptr) {
            subvolume_heat_flux_[static_cast<size_t>(sv)] = phonon_->normalize_to_energy_density(subvolume_heat_flux_[static_cast<size_t>(sv)]);
        }
        subvolume_heat_flux_[static_cast<size_t>(sv)] = mul(subvolume_heat_flux_[static_cast<size_t>(sv)], evpsa2_to_wm2_);
    }

    int axis = 0;
    if (!args_.subvolume_layout.empty() && args_.subvolume_layout.front() == "slice" && args_.subvolume_layout.size() >= 3) {
        axis = std::clamp(std::stoi(args_.subvolume_layout[2]), 0, 2);
    }

    std::vector<double> phi(static_cast<size_t>(nsv), 0.0);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        phi[static_cast<size_t>(sv)] = subvolume_heat_flux_[static_cast<size_t>(sv)][axis];
    }

    const int total_np = std::max(1, std::accumulate(subvolume_particle_counts_.begin(), subvolume_particle_counts_.end(), 0));
    double phi_weighted = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:phi_weighted)
#endif
    for (int sv = 0; sv < nsv; ++sv) {
        phi_weighted += phi[static_cast<size_t>(sv)] * static_cast<double>(subvolume_particle_counts_[static_cast<size_t>(sv)]);
    }
    average_heat_flux_along_axis_ = phi_weighted / static_cast<double>(total_np);

    if (!args_.compute_thermal_conductivity) {
        thermal_conductivity_fit_ = 0.0;
        thermal_conductivity_endpoints_ = 0.0;
        thermal_conductivity_ = 0.0;
        return;
    }

    std::vector<double> T(static_cast<size_t>(nsv + 2), 0.0);
    for (int sv = 0; sv < nsv; ++sv) {
        T[static_cast<size_t>(sv + 1)] = subvolume_temperatures_[static_cast<size_t>(sv)];
    }
    double T_left = subvolume_temperatures_.front();
    double T_right = subvolume_temperatures_.back();
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

    // Method 2: linear fit on middle subvolume_layout.
    const auto& centers = geometry.subvolume_centers();
    int i0 = 0;
    int i1 = nsv - 1;
    if (nsv >= 6) {
        const int cut = std::max(1, nsv / 5);  // remove edge-biased zones near reservoirs
        i0 = cut;
        i1 = nsv - 1 - cut;
    }
    if (i1 - i0 + 1 < 2) {
        i0 = 0;
        i1 = nsv - 1;
    }

    double sx = 0.0;
    double st = 0.0;
    double sxx = 0.0;
    double sxt = 0.0;
    int nfit = 0;
    for (int i = i0; i <= i1 && i < static_cast<int>(centers.size()); ++i) {
        const double x = centers[static_cast<size_t>(i)][axis] * angstrom_to_meter_;
        const double tt = subvolume_temperatures_[static_cast<size_t>(i)];
        sx += x;
        st += tt;
        sxx += x * x;
        sxt += x * tt;
        ++nfit;
    }
    double grad_fit = 0.0;
    const double dT_window = std::abs(
        subvolume_temperatures_[static_cast<size_t>(std::clamp(i1, 0, nsv - 1))] -
        subvolume_temperatures_[static_cast<size_t>(std::clamp(i0, 0, nsv - 1))]);
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
            update_collision_cache(geometry, std::vector<int>{i});
            if (timesteps_to_collision_[i] * time_step_ < 1e-12) {
                particle_positions_[i] = add(particle_positions_[i], mul(particle_velocities_[i], 1e-12));
                update_collision_cache(geometry, std::vector<int>{i});
            }
        }
    }
    (void) phonon;
}

void MonteCarloSolver::write_convergence_header() {
    if (args_.results_base_folder.empty()) {
        return;
    }
    std::ofstream out(std::filesystem::path(args_.results_base_folder) / "convergence.txt", std::ios::trunc);
    out << "# timestep time_ps"; 
    const int nsv = static_cast<int>(subvolume_temperatures_.size());
    for (int i = 0; i < nsv; ++i) {
        out << " T_sv_" << i;
    }
    out << " heatflux kappa_fit kappa_end\n";
}

void MonteCarloSolver::append_convergence_row() const {
    if (args_.results_base_folder.empty()) {
        return;
    }
    std::ofstream out(std::filesystem::path(args_.results_base_folder) / "convergence.txt", std::ios::app);
    // 修改点：在第一列后增加 elapsed_time_
    out << current_timestep_ << " " << elapsed_time_; 
    for (double tsv : subvolume_temperatures_) {
        out << " " << tsv;
    }
    out << " " << average_heat_flux_along_axis_ << " " << thermal_conductivity_fit_ << " " << thermal_conductivity_endpoints_ << '\n';
}

void MonteCarloSolver::run_timestep() {
    if (geometry_ == nullptr) {
        throw std::runtime_error("MonteCarloSolver geometry is not set.");
    }
    if (phonon_ == nullptr) {
        throw std::runtime_error("MonteCarloSolver phonon is not set.");
    }
    const SimulationDomain& geometry = *geometry_;
    const PhononMaterial& phonon = *phonon_;

    const int n_before = particle_count_;
    for (int i = 0; i < n_before; ++i) {
        advance_particle(geometry, phonon, i, time_step_);
        if (i < particle_count_ && i < static_cast<int>(particle_alive_flags_.size()) && particle_alive_flags_[static_cast<size_t>(i)] != 0) {
            particle_subvolume_id_[i] = nearest_subvolume_index(geometry, particle_positions_[i]);
        }
    }
    remove_absorbed_particles();

    const auto injected = inject_particles_from_reservoirs(geometry, phonon);
    if (!injected.empty()) {
        std::vector<int> new_idx;
        new_idx.reserve(injected.size());
        for (const auto& [idx, _] : injected) {
            new_idx.push_back(idx);
        }
        update_collision_cache(geometry, new_idx);
        for (const auto& [idx, dt_in] : injected) {
            const double remain = std::max(0.0, time_step_ - dt_in);
            if (remain > 1e-14 && idx < particle_count_) {
                advance_particle(geometry, phonon, idx, remain);
            }
        }
        remove_absorbed_particles();
    }

    update_particle_temperatures(geometry, phonon);
    apply_lifetime_scattering(phonon);

    ++current_timestep_;
    elapsed_time_ += time_step_;
    if (current_timestep_ % convergence_write_interval_ == 0 || current_timestep_ == 1) {
        update_heat_flux_and_conductivity(geometry);
        append_convergence_row();
    }
    int total_iters = args_.iterations;
    int print_interval = std::max(1, total_iters / 100);

    if (current_timestep_ % print_interval == 0 || current_timestep_ == total_iters) {
        double progress = (static_cast<double>(current_timestep_) / total_iters) * 100.0;
        
        std::cout << "--- Progress: " << std::fixed << std::setprecision(1) << progress << "% ---" << std::endl;
        
        std::cout << "Temperature Profile (K): ";
        for (double t : subvolume_temperatures_) {
            std::cout << std::setprecision(2) << t << " ";
        }
        std::cout << std::endl;

        if (args_.compute_thermal_conductivity) {
            std::cout << "Current Conductivity (Fit): " << thermal_conductivity_fit_ << " W/mK" << std::endl;
            std::cout << "Current Conductivity (End): " << thermal_conductivity_endpoints_ << " W/mK" << std::endl;
        }
        std::cout << "-------------------------" << std::endl;
    }
}
