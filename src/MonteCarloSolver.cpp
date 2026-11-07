#include "MonteCarloSolver.h"

#include "PhononMaterial.h"
#include "SimulationDomain.h"
#include "material/ScatteringMatrixReader.h"
#include "solver/StageTimer.h"

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

// 函数说明：执行三维向量减法，服务于碰撞几何与距离计算。
MonteCarloSolver::Vec3 MonteCarloSolver::sub(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

// 函数说明：计算向量点积，用于投影、入射角与法向分量判断。
double MonteCarloSolver::dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// 函数说明：提供线程独立随机数发生器，保证并行采样过程的线程安全。
std::mt19937_64& MonteCarloSolver::thread_rng() const {
#ifdef PHONOMC_USE_OPENMP
    return thread_rngs_[static_cast<size_t>(omp_get_thread_num())];
#else
    return rng_;
#endif
}

// 函数说明：构造蒙特卡洛求解器并完成粒子、边界、热流统计与输出初始化。
MonteCarloSolver::MonteCarloSolver(const SimulationConfig& args,
                                   const SimulationDomain& geometry,
                                   const PhononMaterial& phonon)
    : MonteCarloSolver(
          args,
          geometry,
          std::vector<std::reference_wrapper<const PhononMaterial>>{std::cref(phonon)}) {}

MonteCarloSolver::MonteCarloSolver(
    const SimulationConfig& args,
    const SimulationDomain& geometry,
    const std::vector<std::reference_wrapper<const PhononMaterial>>& materials)
    : args_(args), result_writer_(args.output_folder), convergence_monitor_(args),
      geometry_(&geometry) {
    if (materials.empty()) {
        throw std::runtime_error("MonteCarloSolver requires at least one material.");
    }
    materials_.reserve(materials.size());
    for (const auto& entry : materials) {
        materials_.push_back(&entry.get());
    }
    phonon_ = materials_.front();
    validate_simulation_controls(args_);
    validate_scattering_config(args_);
    if (!args_.material_folders.empty() && args_.material_folders.size() != materials_.size()) {
        throw std::runtime_error("Loaded material count does not match materials.folders.");
    }
    using Clock = std::chrono::steady_clock;
    const auto t_init_begin = Clock::now();
    rng_.seed(args_.random_seed);
    initial_particle_count_ = static_cast<int>(std::llround(args_.particle_count));
    particles_.reset(initial_particle_count_);
    push_eps_ = 1e-10 * std::max(args_.time_step, 1.0);
    particle_density_ = static_cast<double>(particles_.size()) / geometry.volume();
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
        for (const auto* m : materials_) {
            densities.push_back(m->active_mode_count() / m->energy_density_normalization());
        }
        double states = 0.0;
        for (size_t cell = 0; cell < grid_.material_ids.size(); ++cell) {
            states += geometry.grid_volumes()[cell] * densities[grid_.material_ids[cell]];
        }
        const double common_weight = states / initial_particle_count_;
        for (double density : densities) {
            if (!(density > 0.0) || !std::isfinite(density)) {
                throw std::runtime_error("Invalid material phase-space density.");
            }
            particles_.material_spatial_volumes_a3.push_back(common_weight / density);
        }
    }
    if (args_.collision_model == CollisionModel::FullMatrix) {
        std::vector<phonomc::ScatteringMatrix> matrices;
        for (const auto& name : args_.scattering_matrix_files) {
            std::filesystem::path file(name);
            if (file.is_relative()) {
                file = std::filesystem::path(args_.input_directory) / file;
            }
            matrices.emplace_back(phonomc::read_scattering_matrix(file));
        }
        full_matrix_collision_.configure(std::move(matrices), materials_, background_temperature_);
    } else if (args_.collision_model == CollisionModel::Callaway) {
        if (args_.callaway_nonlinear) {
            nonlinear_callaway_.configure(materials_, background_temperature_);
        } else {
            carrier_collision_.configure(materials_, background_temperature_);
        }
    } else if (args_.temperature_gradient) {
        carrier_collision_.configure(materials_, background_temperature_, true);
    }
    background_cache_.reserve(materials_.size());
    for (size_t mid = 0; mid < materials_.size(); ++mid) {
        background_cache_.emplace_back(
            *materials_[mid],
            background_temperature_,
            particles_.spatial_volume(static_cast<int>(mid), particle_spatial_weight_a3_),
            uses_linearized_transport(args_));
    }

    std::cout << "MonteCarloSolver initialized: particle_count=" << particles_.size()
              << ", time_step=" << args_.time_step
              << ", convergence_write_interval=" << args_.convergence_write_interval
              << ", progress_temperature_summary_only="
              << (args_.progress_temperature_summary_only ? "true" : "false")
              << ", density=" << particle_density_ << '\n';
    std::cout << "Thermal conductivity estimation: "
              << (args_.compute_kappa ? "enabled" : "disabled") << '\n';
    std::cout << "Temperature reference: background=fixed T=" << background_temperature_ << " K";
    std::cout << ", collision="
              << (args_.collision_model == CollisionModel::Callaway
                      ? (args_.callaway_nonlinear
                             ? "Callaway N/R (nonlinear displaced Bose)"
                             : "Callaway N/R (U + supplied disorder; fixed reference)")
                  : args_.collision_model == CollisionModel::FullMatrix
                      ? "full_matrix (fixed reference)"
                      : "RTA (existing carriers)");
    if (args_.collision_model == CollisionModel::Rta) {
        std::cout << ", lifetime=" << (fixed_lifetime_temperature_ ? "fixed" : "local");
    }
    if (args_.collision_model == CollisionModel::Rta && fixed_lifetime_temperature_) {
        std::cout << " T=" << lifetime_temperature_ << " K";
    }
    std::cout << '\n';
#ifdef PHONOMC_USE_OPENMP
    openmp_thread_count_ = std::max(1, omp_get_max_threads());
    thread_rngs_.resize(static_cast<size_t>(openmp_thread_count_));
    for (int tid = 0; tid < openmp_thread_count_; ++tid) {
        std::seed_seq sequence{static_cast<unsigned>(args_.random_seed & 0xffffffffu),
                               static_cast<unsigned>((args_.random_seed >> 32) & 0xffffffffu),
                               static_cast<unsigned>(tid),
                               0x504d4301u};
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
        if (materials_.size() != 1) {
            throw std::runtime_error("Gradient driving requires one material.");
        }
        gradient_drive_.configure(args_, geometry, *phonon_, particle_spatial_weight_a3_);
        std::cout << "Periodic linear response: gradient=" << *args_.temperature_gradient
                  << " K/m, transport_axis=" << "xyz"[args_.transport_axis]
                  << "; conductivity=-mean_flux/imposed_gradient; T output is T0 + periodic "
                     "correction.\n";
    }
    initialize_interface_mode_banks();
    initialize_particles(geometry);
    initialize_local_heat_source(geometry);
    if (args_.gradient_warm_start_ps > 0) {
        grid_.cells.ensure(particles_, static_cast<int>(grid_.material_ids.size()));
        auto seeded =
            gradient_drive_.warm_start(particles_, grid_.cells, args_.gradient_warm_start_ps, rng_);
        grid_.cells.invalidate();
        transport_->update_collision_cache(particles_, seeded.appended_indices);
        refresh_thermodynamic_state(geometry);
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
    return material_layout_.material_at(position);
}

double MonteCarloSolver::material_energy_weight(int material_id) const {
    return background_cache_.at(static_cast<size_t>(material_id)).energy_weight();
}

std::vector<double> MonteCarloSolver::material_energy_weights() const {
    std::vector<double> weights;
    weights.reserve(materials_.size());
    for (size_t material_id = 0; material_id < materials_.size(); ++material_id) {
        weights.push_back(material_energy_weight(static_cast<int>(material_id)));
    }
    return weights;
}

void MonteCarloSolver::apply_phonon_collisions() {
    phonomc::StageTimer timer(timer_lifetime_, profile_timers_enabled_);
    if (args_.callaway_nonlinear) {
        energy_ledger_.set_lifetime_residual(nonlinear_callaway_.apply(
            particles_,
            materials_,
            {grid_.cells.counts(), grid_.material_ids, grid_.temperatures, &grid_.cells},
            {args_.time_step,
             particle_spatial_weight_a3_,
             !fixed_lifetime_temperature_,
             lifetime_temperature_}));
        return;
    }
    if (args_.collision_model == CollisionModel::Callaway || args_.temperature_gradient) {
        energy_ledger_.set_lifetime_residual(carrier_collision_.apply(
            particles_,
            materials_,
            {grid_.cells.counts(), grid_.material_ids, grid_.temperatures, &grid_.cells},
            args_.time_step,
            particle_spatial_weight_a3_));
        return;
    }
    if (args_.collision_model == CollisionModel::FullMatrix) {
        const auto result = full_matrix_collision_.apply(
            particles_,
            materials_,
            {grid_.cells.counts(), grid_.material_ids, grid_.temperatures, &grid_.cells},
            args_.time_step,
            particle_spatial_weight_a3_,
            background_temperature_,
            rng_,
            args_.temperature_gradient.has_value());
        energy_ledger_.set_lifetime_residual(result.residual_ev);
        if (!result.appended_indices.empty()) {
            grid_.cells.invalidate();
            grid_.cells.ensure(particles_, static_cast<int>(grid_.material_ids.size()));
            transport_->update_collision_cache(particles_, result.appended_indices);
        }
        return;
    }
    energy_ledger_.set_lifetime_residual(rta_collision_.apply(
        particles_,
        materials_,
        {grid_.cells.counts(), grid_.material_ids, grid_.temperatures, &grid_.cells},
        {args_.time_step,
         particle_spatial_weight_a3_,
         !fixed_lifetime_temperature_,
         lifetime_temperature_}));
}

void MonteCarloSolver::apply_temperature_gradient(double time_ps) {
    if (!args_.temperature_gradient) {
        return;
    }
    grid_.cells.ensure(particles_, static_cast<int>(grid_.material_ids.size()));
    const auto result = gradient_drive_.apply(particles_, grid_.cells, time_ps, rng_);
    energy_ledger_.record_drive(result.occupation_updates, result.energy_ev);
    if (!result.appended_indices.empty()) {
        grid_.cells.invalidate();
        grid_.cells.ensure(particles_, static_cast<int>(grid_.material_ids.size()));
        transport_->update_collision_cache(particles_, result.appended_indices);
    }
}

void MonteCarloSolver::transport_particles() {
    // All transport/insertion/compaction changes occur before the next refresh.
    grid_.cells.invalidate();
    {
        phonomc::StageTimer timer(timer_advance_main_, profile_timers_enabled_);
        transport_->advance_existing(particles_, boundary_ledger_, {rng_, thread_rngs_});
    }
    {
        phonomc::StageTimer timer(timer_remove_absorb_1_, profile_timers_enabled_);
        (void)particles_.compact_alive();
    }

    std::vector<std::pair<int, double>> injected;
    {
        phonomc::StageTimer timer(timer_inject_build_, profile_timers_enabled_);
        injected = reservoir_injector_->inject(particles_, boundary_ledger_, thread_rng());
    }
    if (!injected.empty()) {
        {
            phonomc::StageTimer timer(timer_inject_cache_, profile_timers_enabled_);
            std::vector<int> new_idx;
            new_idx.reserve(injected.size());
            for (const auto& [idx, _] : injected) {
                new_idx.push_back(idx);
            }
            transport_->update_collision_cache(particles_, new_idx);
        }
        {
            phonomc::StageTimer timer(timer_advance_injected_, profile_timers_enabled_);
            transport_->advance_injected(particles_, injected, boundary_ledger_, {rng_, thread_rngs_});
        }
        {
            phonomc::StageTimer timer(timer_remove_absorb_2_, profile_timers_enabled_);
            (void)particles_.compact_alive();
        }
    }

    {
        phonomc::StageTimer timer(timer_boundary_checks_, profile_timers_enabled_);
        boundary_ledger_.finish_transport();
        if (escaped_recovery_check_interval_ > 0 &&
            (current_timestep_ % escaped_recovery_check_interval_) == 0) {
            transport_->validate_containment(particles_);
        }
    }
}

void MonteCarloSolver::run_timestep() {
    if (incomplete_timestep_) {
        throw std::runtime_error(
            "Cannot resume a failed timestep; construct a new solver after fixing its cause.");
    }
    assert(particles_.aligned());
    incomplete_timestep_ = true;
    {
        phonomc::StageTimer total(timer_total_, profile_timers_enabled_);
        energy_ledger_.begin_step();
        transport_->begin_step();
#ifdef PHONOMC_USE_OPENMP
        boundary_ledger_.begin_step(omp_get_max_threads());
#else
        boundary_ledger_.begin_step(1);
#endif
        const double midpoint = elapsed_time_ + 0.5 * args_.time_step;
        const double end = elapsed_time_ + args_.time_step;

        // Symmetric source splitting: half the drive on each side of transport
        // and collision. Heating must refresh the local scattering temperature.
        apply_heat_source(elapsed_time_, midpoint);
        if (energy_ledger_.step().source_energy_ev > 0.0) {
            refresh_thermodynamic_state(*geometry_);
        }
        apply_temperature_gradient(0.5 * args_.time_step);

        transport_particles();
        refresh_thermodynamic_state(*geometry_);
        apply_phonon_collisions();

        apply_heat_source(midpoint, end);
        apply_temperature_gradient(0.5 * args_.time_step);
        refresh_thermodynamic_state(*geometry_);

        // Close the energy balance only after the final temperature reconstruction.
        {
            phonomc::StageTimer accounting(timer_energy_ledger_, profile_timers_enabled_);
            energy_ledger_.finish_step(compute_total_thermal_energy_ev(*geometry_),
                                       boundary_ledger_.step());
            boundary_ledger_.commit_step();
        }
        ++current_timestep_;
        elapsed_time_ += args_.time_step;
        record_observables_if_due();
    }
    report_progress_if_needed();
    report_timestep_timers_if_needed();
    assert(particles_.aligned());
    incomplete_timestep_ = false;
}
