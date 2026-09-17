#include "MonteCarloSolver.h"

#include "PhononMaterial.h"
#include "SimulationDomain.h"
#include "material/ScatteringMatrixReader.h"
#include "solver/ModalResampler.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

#ifdef PHONOMC_USE_OPENMP
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

// 函数说明：提供线程独立随机数发生器，保证并行采样过程的线程安全。
std::mt19937_64& MonteCarloSolver::thread_rng() const {
#ifdef PHONOMC_USE_OPENMP
    const int tid = omp_get_thread_num();
    if (tid >= 0 && tid < static_cast<int>(thread_rngs_.size())) {
        return thread_rngs_[static_cast<size_t>(tid)];
    }
    return rng_;
#else
    return rng_;
#endif
}

// 函数说明：构造蒙特卡洛求解器并完成粒子、边界、热流统计与输出初始化。
MonteCarloSolver::MonteCarloSolver(
    const SimulationConfig& args,
    const SimulationDomain& geometry,
    const PhononMaterial& phonon)
    : MonteCarloSolver(args, geometry, std::vector<std::reference_wrapper<const PhononMaterial>> {std::cref(phonon)}) {}

MonteCarloSolver::MonteCarloSolver(
    const SimulationConfig& args,
    const SimulationDomain& geometry,
    const std::vector<std::reference_wrapper<const PhononMaterial>>& materials)
    : args_(args), result_writer_(args.output_folder), convergence_monitor_(args), geometry_(&geometry) {
    if (materials.empty()) {
        throw std::runtime_error("MonteCarloSolver requires at least one material.");
    }
    materials_.reserve(materials.size());
    for (const auto& entry : materials) {
        materials_.push_back(&entry.get());
    }
    phonon_ = materials_.front();
    validate_scattering_config(args_);
    if (!args_.material_folders.empty() && args_.material_folders.size() != materials_.size()) {
        throw std::runtime_error("Loaded material count does not match materials.folders.");
    }
    using Clock = std::chrono::steady_clock;
    const auto t_init_begin = Clock::now();
    rng_seed_base_ = args_.random_seed;
    rng_.seed(rng_seed_base_);
    resampling_rng_.seed(rng_seed_base_ ^ 0x524553414d504c45ULL);
    initial_particle_count_ = std::max(1, static_cast<int>(std::llround(args_.particle_count)));
    particles_.reset(initial_particle_count_);
    time_step_ = std::max(1e-12, args_.time_step);
    convergence_write_interval_ = std::max(1, args_.convergence_write_interval);
    push_eps_ = 1e-10 * std::max(time_step_, 1.0);
    particle_density_ = static_cast<double>(particles_.size()) / std::max(geometry.volume(), 1e-12);
    particle_spatial_weight_a3_ = geometry.volume() / static_cast<double>(initial_particle_count_);
    fixed_lifetime_temperature_ = !args_.lifetime_temperature_is_local;
    background_temperature_ = args_.background_temperature;
    lifetime_temperature_ = args_.lifetime_temperature;
    initialize_material_layout(geometry);
    if (materials_.size() > 1) {
        // Uniform carrier density gives different state weights in different
        // crystals. DMM mixes physical phase-space fluxes, so sample carriers
        // in proportion to the active state density and use a common weight.
        std::vector<double> densities;
        for (const auto* m : materials_)
            densities.push_back(m->active_mode_count()/m->energy_density_normalization());
        double states = 0.0;
        for (size_t cell=0; cell<grid_.material_ids.size(); ++cell)
            states += geometry.grid_volumes()[cell]*densities[grid_.material_ids[cell]];
        const double common_weight = states/initial_particle_count_;
        for (double density : densities) {
            if (!(density > 0.0) || !std::isfinite(density))
                throw std::runtime_error("Invalid material phase-space density.");
            particles_.material_spatial_volumes_a3.push_back(common_weight/density);
        }
    }
    if (args_.collision_model == CollisionModel::FullMatrix) {
        std::vector<phonomc::ScatteringMatrix> matrices;
        for (const auto& name : args_.scattering_matrix_files) {
            std::filesystem::path file(name);
            if (file.is_relative()) file = std::filesystem::path(args_.input_directory) / file;
            matrices.emplace_back(phonomc::read_scattering_matrix(file));
        }
        full_matrix_collision_.configure(std::move(matrices), materials_, background_temperature_);
    } else if (args_.collision_model == CollisionModel::Callaway) {
        if(args_.callaway_nonlinear) nonlinear_callaway_.configure(materials_,background_temperature_);
        else carrier_collision_.configure(materials_,background_temperature_);
    } else if (args_.temperature_gradient) {
        carrier_collision_.configure(materials_,background_temperature_,true);
    }
    background_cache_.reserve(materials_.size());
    for (size_t mid=0; mid<materials_.size(); ++mid)
        background_cache_.emplace_back(*materials_[mid], background_temperature_,
                                       particles_.spatial_volume(static_cast<int>(mid),particle_spatial_weight_a3_),
                                       uses_linearized_transport(args_));

    std::cout << "MonteCarloSolver initialized: particle_count=" << particles_.size()
              << ", time_step=" << time_step_
              << ", convergence_write_interval=" << convergence_write_interval_
              << ", progress_temperature_summary_only="
              << (args_.progress_temperature_summary_only ? "true" : "false")
              << ", density=" << particle_density_ << '\n';
    std::cout << "Thermal conductivity estimation: "
              << (args_.compute_kappa ? "enabled" : "disabled")
              << '\n';
    std::cout << "Temperature reference: background=fixed T="
              << background_temperature_ << " K";
    std::cout << ", collision=" << (args_.collision_model == CollisionModel::Callaway ? (args_.callaway_nonlinear ? "Callaway N/R (nonlinear displaced Bose)" : "Callaway N/R (U + supplied disorder; fixed reference)") :
        args_.collision_model == CollisionModel::FullMatrix ? "full_matrix (fixed reference)" :
        "RTA (existing carriers)");
    if (args_.collision_model == CollisionModel::Rta) std::cout << ", lifetime=" << (fixed_lifetime_temperature_ ? "fixed" : "local");
    if (args_.collision_model == CollisionModel::Rta && fixed_lifetime_temperature_) {
        std::cout << " T=" << lifetime_temperature_ << " K";
    }
    std::cout << '\n';
#ifdef PHONOMC_USE_OPENMP
    openmp_thread_count_ = std::max(1, omp_get_max_threads());
    thread_rngs_.resize(static_cast<size_t>(openmp_thread_count_));
    for (int tid = 0; tid < openmp_thread_count_; ++tid) {
        std::seed_seq sequence {
            static_cast<unsigned>(rng_seed_base_ & 0xffffffffu),
            static_cast<unsigned>((rng_seed_base_ >> 32) & 0xffffffffu),
            static_cast<unsigned>(tid),
            0x504d4301u
        };
        thread_rngs_[static_cast<size_t>(tid)].seed(sequence);
    }
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

    if (args_.temperature_gradient) {
        if (materials_.size()!=1) throw std::runtime_error("Gradient driving requires one material.");
        gradient_drive_.configure(args_,geometry,*phonon_,particle_spatial_weight_a3_);
        std::cout << "Periodic linear response: gradient=" << *args_.temperature_gradient
                  << " K/m, transport_axis=" << "xyz"[args_.transport_axis]
                  << "; conductivity=-mean_flux/imposed_gradient; T output is T0 + periodic correction.\n";
    }
    initialize_interface_mode_banks();
    initialize_particles(geometry);
    initialize_local_heat_source(geometry);
    if (args_.gradient_warm_start_ps>0) {
        grid_.cells.ensure(particles_,static_cast<int>(grid_.material_ids.size()));
        auto seeded=gradient_drive_.warm_start(particles_,grid_.cells,args_.gradient_warm_start_ps,rng_);
        grid_.cells.invalidate();
        update_collision_cache(geometry,seeded.appended_indices);
        update_particle_temperatures(geometry);
    }
    energy_ledger_.initialize(compute_total_thermal_energy_ev(geometry));
    write_convergence_header();
    update_heat_flux_and_conductivity(geometry);
    sample_convergence();
    append_convergence_row();
    const auto t_init_end = Clock::now();
    const double init_sec = std::chrono::duration<double>(t_init_end - t_init_begin).count();
    std::cout << "[init] Initialization complete in " << std::fixed << std::setprecision(2)
              << init_sec << " s\n";
}

const PhononMaterial& MonteCarloSolver::material(int material_id) const {
    if (material_id < 0 || material_id >= static_cast<int>(materials_.size())) {
        throw std::runtime_error("Particle references an invalid material index.");
    }
    return *materials_[static_cast<size_t>(material_id)];
}

const PhononMaterial& MonteCarloSolver::particle_material(int particle_index) const {
    return material(particles_.material_ids.at(static_cast<size_t>(particle_index)));
}

int MonteCarloSolver::material_index_at(const Vec3& position) const {
    if (materials_.size() <= 1 || args_.material_interface_positions.empty()) {
        return 0;
    }
    const double coordinate = position[static_cast<size_t>(args_.material_interface_axis)];
    return static_cast<int>(std::upper_bound(
        args_.material_interface_positions.begin(),
        args_.material_interface_positions.end(),
        coordinate) - args_.material_interface_positions.begin());
}

double MonteCarloSolver::material_energy_weight(int material_id) const {
    return background_cache_.at(static_cast<size_t>(material_id)).energy_weight();
}

// 函数说明：返回全局不变的偏差能量背景温度。
double MonteCarloSolver::background_temperature_reference() const {
    return background_temperature_;
}

void MonteCarloSolver::apply_lifetime_scattering() {
    if(args_.callaway_nonlinear) {
        energy_ledger_.set_lifetime_residual(nonlinear_callaway_.apply(particles_,materials_,
            {grid_.cells.counts(),grid_.material_ids,grid_.temperatures,&grid_.cells},
            {time_step_,particle_spatial_weight_a3_,!fixed_lifetime_temperature_,lifetime_temperature_}));
        return;
    }
    if (args_.collision_model == CollisionModel::Callaway || args_.temperature_gradient) {
        energy_ledger_.set_lifetime_residual(carrier_collision_.apply(particles_,materials_,
            {grid_.cells.counts(),grid_.material_ids,grid_.temperatures,&grid_.cells},
            time_step_,particle_spatial_weight_a3_));
        return;
    }
    if (args_.collision_model == CollisionModel::FullMatrix) {
        const auto result = full_matrix_collision_.apply(particles_, materials_,
            {grid_.cells.counts(), grid_.material_ids, grid_.temperatures, &grid_.cells},
            time_step_, particle_spatial_weight_a3_, background_temperature_, rng_, args_.temperature_gradient.has_value());
        energy_ledger_.set_lifetime_residual(result.residual_ev);
        if (!result.appended_indices.empty()) {
            grid_.cells.invalidate();
            grid_.cells.ensure(particles_, static_cast<int>(grid_.material_ids.size()));
            update_collision_cache(*geometry_, result.appended_indices);
        }
        return;
    }
    energy_ledger_.set_lifetime_residual(rta_collision_.apply(
        particles_, materials_,
        {grid_.cells.counts(), grid_.material_ids, grid_.temperatures, &grid_.cells},
        {time_step_, particle_spatial_weight_a3_, !fixed_lifetime_temperature_, lifetime_temperature_}));
}

void MonteCarloSolver::apply_temperature_gradient(double time_ps) {
    if (!args_.temperature_gradient) return;
    grid_.cells.ensure(particles_,static_cast<int>(grid_.material_ids.size()));
    const auto result=gradient_drive_.apply(particles_,grid_.cells,time_ps,rng_);
    energy_ledger_.record_drive(result.occupation_updates,result.energy_ev);
    if (!result.appended_indices.empty()) {
        grid_.cells.invalidate();
        grid_.cells.ensure(particles_,static_cast<int>(grid_.material_ids.size()));
        update_collision_cache(*geometry_,result.appended_indices);
    }
}

// 函数说明：执行一次完整时间步流程：推进、注入、温度更新、散射与输出。
void MonteCarloSolver::run_timestep() {
    if (incomplete_timestep_)
        throw std::runtime_error("Cannot resume a failed timestep; construct a new solver after fixing its cause.");
    if (geometry_ == nullptr) {
        throw std::runtime_error("MonteCarloSolver geometry is not set.");
    }
    if (materials_.empty()) {
        throw std::runtime_error("MonteCarloSolver materials are not set.");
    }
    assert(particles_.aligned());
    incomplete_timestep_ = true;
    const SimulationDomain& geometry = *geometry_;
    using Clock = std::chrono::steady_clock;
    const auto t_step_begin = Clock::now();
    energy_ledger_.begin_step();
    excessive_collision_particle_.store(-1, std::memory_order_relaxed);
#ifdef PHONOMC_USE_OPENMP
    boundary_ledger_.begin_step(std::max(1, omp_get_max_threads()));
#else
    boundary_ledger_.begin_step(1);
#endif

    const double midpoint_time = elapsed_time_ + 0.5 * time_step_;
    const double first_source_time_ps = local_heat_source_integrated_time_factor(
        elapsed_time_, midpoint_time);
    const double second_source_time_ps = local_heat_source_integrated_time_factor(
        midpoint_time, elapsed_time_ + time_step_);
    if (heat_source_.enabled() && (first_source_time_ps > 0.0 || heat_source_.pending_energy_ev() > 0.0)) {
        apply_local_heat_source_to_occupations(first_source_time_ps);
        // First half of a Strang-split source must update the local reference
        // temperature before transport and boundary scattering.
        if (energy_ledger_.step().source_energy_ev > 0.0) update_particle_temperatures(geometry);
    }

    apply_temperature_gradient(0.5*time_step_);
    // All transport/insertion/compaction changes occur before the next refresh.
    grid_.cells.invalidate();
    const int n_before = particles_.size();
    const auto t_advance_begin = Clock::now();
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (int i = 0; i < n_before; ++i) {
        advance_particle(geometry, i, time_step_);
        if (i < particles_.size() && i < static_cast<int>(particles_.alive.size()) && particles_.alive[static_cast<size_t>(i)] != 0) {
            particles_.grid_ids[i] = nearest_grid_index(geometry, particles_.positions[i]);
        }
    }
    throw_if_excessive_collisions();
    throw_if_collision_cache_failed("advancing particles");
    const auto t_advance_end = Clock::now();
    const auto t_remove1_begin = t_advance_end;
    (void) particles_.compact_alive();
    const auto t_remove1_end = Clock::now();

    const auto t_inject_begin = t_remove1_end;
    const auto injected_reservoir = inject_particles_from_reservoirs(geometry);
    const auto& injected = injected_reservoir;
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
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (int k = 0; k < ninj; ++k) {
            const auto [idx, dt_in] = injected[static_cast<size_t>(k)];
            const double remain = std::max(0.0, time_step_ - dt_in);
            if (remain > 1e-14 && idx < particles_.size()) {
                advance_particle(geometry, idx, remain);
            }
            if (idx < particles_.size() &&
                idx < static_cast<int>(particles_.alive.size()) &&
                particles_.alive[static_cast<size_t>(idx)] != 0) {
                particles_.grid_ids[idx] = nearest_grid_index(geometry, particles_.positions[idx]);
            }
        }
        throw_if_excessive_collisions();
        throw_if_collision_cache_failed("advancing injected particles");
        const auto t_adv_inj_end = Clock::now();
        const auto t_remove2_begin = t_adv_inj_end;
        (void) particles_.compact_alive();
        const auto t_remove2_end = Clock::now();
        if (profile_timers_enabled_) {
            timer_inject_cache_ += std::chrono::duration<double>(t_cache_end - t_cache_begin).count();
            timer_advance_injected_ += std::chrono::duration<double>(t_adv_inj_end - t_adv_inj_begin).count();
            timer_remove_absorb_2_ += std::chrono::duration<double>(t_remove2_end - t_remove2_begin).count();
        }
    }

    const auto t_boundary_begin = Clock::now();
    boundary_ledger_.finish_transport();

    if (escaped_recovery_check_interval_ > 0 &&
        (current_timestep_ % escaped_recovery_check_interval_) == 0) {
        throw_if_particles_escaped(geometry);
    }

    const auto t_boundary_end = Clock::now();
    update_particle_temperatures(geometry);
    const auto t_life_begin = Clock::now();
    apply_lifetime_scattering();
    const auto t_life_end = Clock::now();

    if (heat_source_.enabled() && (second_source_time_ps > 0.0 || heat_source_.pending_energy_ev() > 0.0)) {
        apply_local_heat_source_to_occupations(second_source_time_ps);
    }
    apply_temperature_gradient(0.5*time_step_);
    // Occupations changed during collision and the second source half. Refresh
    // particle energies and temperatures before heat-flux/output statistics.
    update_particle_temperatures(geometry);
    const auto t_ledger_begin = Clock::now();
    if (args_.resample_interval>0 && (current_timestep_+1)%args_.resample_interval==0) {
        const auto start=Clock::now();
        const auto res=phonomc::resample_modal_carriers(particles_,grid_.cells,
            static_cast<int>(grid_.material_ids.size()),*phonon_,background_cache_.front(),
            args_.resample_per_mode_sign,resampling_rng_);
        resampling_removed_+=res.removed;resampling_residual_ev_+=res.energy_residual_ev;
        update_particle_temperatures(geometry);
        resampling_seconds_+=std::chrono::duration<double>(Clock::now()-start).count();
    }
    energy_ledger_.finish_step(compute_total_thermal_energy_ev(geometry), boundary_ledger_.step());
    boundary_ledger_.commit_step();

    ++current_timestep_;
    elapsed_time_ += time_step_;
    const auto t_stats_begin = Clock::now();
    if (current_timestep_ % convergence_write_interval_ == 0 || current_timestep_ == 1) {
        update_heat_flux_and_conductivity(geometry);
        sample_convergence();
        append_convergence_row();
    }
    const auto t_stats_end = Clock::now();
    const auto t_step_end = t_stats_end;
    if (profile_timers_enabled_) {
        timer_total_ += std::chrono::duration<double>(t_step_end - t_step_begin).count();
        timer_advance_main_ += std::chrono::duration<double>(t_advance_end - t_advance_begin).count();
        timer_remove_absorb_1_ += std::chrono::duration<double>(t_remove1_end - t_remove1_begin).count();
        timer_inject_build_ += std::chrono::duration<double>(t_inject_end - t_inject_begin).count();
        timer_boundary_checks_ += std::chrono::duration<double>(t_boundary_end - t_boundary_begin).count();
        timer_energy_ledger_ += std::chrono::duration<double>(t_stats_begin - t_ledger_begin).count();
        timer_lifetime_ += std::chrono::duration<double>(t_life_end - t_life_begin).count();
        timer_stats_ += std::chrono::duration<double>(t_stats_end - t_stats_begin).count();
    }
    int total_iters = args_.iterations;
    int print_interval = std::max(1, total_iters / 100);

    if (current_timestep_ % print_interval == 0 || current_timestep_ == total_iters) {
        double progress = (static_cast<double>(current_timestep_) / total_iters) * 100.0;
        
        std::cout << "--- Progress: " << std::fixed << std::setprecision(1) << progress << "% ---" << std::endl;
        
        if (args_.progress_temperature_summary_only) {
            if (grid_.temperatures.empty()) {
                std::cout << "Temperature Summary (K): Tmin=nan, Tavg=nan, Tmax=nan" << std::endl;
            } else {
                const auto mm = std::minmax_element(grid_.temperatures.begin(), grid_.temperatures.end());
                const double tsum = std::accumulate(grid_.temperatures.begin(), grid_.temperatures.end(), 0.0);
                const double tavg = tsum / static_cast<double>(grid_.temperatures.size());
                std::cout << "Temperature Summary (K): Tmin=" << std::setprecision(2) << *mm.first
                          << ", Tavg=" << tavg
                          << ", Tmax=" << *mm.second << std::endl;
            }
        } else {
            std::cout << "Temperature Profile (K): ";
            for (double t : grid_.temperatures) {
                std::cout << std::setprecision(2) << t << " ";
            }
            std::cout << std::endl;
        }

        if (args_.compute_kappa) {
            if (!args_.temperature_gradient)
                std::cout << "Current Conductivity (Int): " << thermal_conductivity_fit_ << " W/mK" << std::endl;
            std::cout << (args_.temperature_gradient ? "Current Conductivity (Gradient): " : "Current Conductivity (Eff): ")
                      << thermal_conductivity_endpoints_ << " W/mK" << std::endl;
        }
        std::cout << "Reservoir Balance (step): absorbed=" << boundary_ledger_.step().absorbed
                  << ", injected=" << boundary_ledger_.step().injected
                  << ", recovered=" << 0
                  << ", net=" << boundary_ledger_.step().net()
                  << ", particle_count=" << particles_.size() << std::endl;
        if (heat_source_.enabled()) {
            std::cout << "Local Heat Source (step): occupation_updates=" << energy_ledger_.step().source_updates
                      << ", injected_energy_ev=" << energy_ledger_.step().source_energy_ev
                      << ", prescribed_total_ev=" << heat_source_.prescribed_energy_ev()
                      << ", pending_total_ev=" << heat_source_.pending_energy_ev() << std::endl;
        }
        std::cout << "Energy Diagnostics (eV): reservoir_absorbed=" << boundary_ledger_.step().absorbed_energy_ev
                  << ", reservoir_injected=" << boundary_ledger_.step().injected_energy_ev
                  << ", lifetime_residual=" << energy_ledger_.step().lifetime_residual_ev
                  << ", balance_residual=" << energy_ledger_.step().balance_residual_ev
                  << ", total_thermal=" << energy_ledger_.thermal_energy_ev() << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
    assert(particles_.aligned());
    report_timestep_timers_if_needed();
    incomplete_timestep_ = false;
}
