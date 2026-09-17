#include "MonteCarloSolver.h"

#include "PhononMaterial.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

void MonteCarloSolver::initialize_interface_mode_banks() {
    if (!std::isfinite(args_.material_interface_frequency_bin_thz) || args_.material_interface_frequency_bin_thz<=0)
        throw std::runtime_error("Interface frequency bin must be finite and positive.");
    interface_sampler_.configure(
        materials_.size() > 1 && !args_.material_interface_positions.empty()
            ? materials_ : std::vector<const PhononMaterial*> {},
        args_.material_interface_axis);
    if (!std::isfinite(args_.material_interface_amm_fraction)||args_.material_interface_amm_fraction<0||args_.material_interface_amm_fraction>1)
        throw std::runtime_error("AMM fraction must be in [0,1].");
    if (args_.material_interface_model==InterfaceModel::AcousticMismatch ||
        (args_.material_interface_model==InterfaceModel::MixedMismatch && args_.material_interface_amm_fraction>0))
        acoustic_interface_sampler_.configure(materials_,args_.material_interface_axis,
            args_.material_mass_densities_kg_m3,args_.material_interface_frequency_bin_thz,
            args_.material_interface_parallel_bin_inv_a);
}

double MonteCarloSolver::next_interface_time(int particle_index, int& next_material) const {
    next_material = -1;
    if (materials_.size() <= 1) {
        return std::numeric_limits<double>::infinity();
    }
    const int current = particles_.material_ids[static_cast<size_t>(particle_index)];
    const int axis = args_.material_interface_axis;
    const double velocity = particles_.velocities[static_cast<size_t>(particle_index)][static_cast<size_t>(axis)];
    if (std::abs(velocity) <= 1e-18) {
        return std::numeric_limits<double>::infinity();
    }
    double plane = 0.0;
    if (velocity > 0.0 && current + 1 < static_cast<int>(materials_.size())) {
        plane = args_.material_interface_positions[static_cast<size_t>(current)];
        next_material = current + 1;
    } else if (velocity < 0.0 && current > 0) {
        plane = args_.material_interface_positions[static_cast<size_t>(current - 1)];
        next_material = current - 1;
    } else {
        return std::numeric_limits<double>::infinity();
    }
    const double time = (plane - particles_.positions[static_cast<size_t>(particle_index)][static_cast<size_t>(axis)]) / velocity;
    if (!std::isfinite(time) || time <= 1e-14) {
        next_material = -1;
        return std::numeric_limits<double>::infinity();
    }
    return time;
}

void MonteCarloSolver::process_material_interface(int particle_index, int next_material) {
    const int old_material = particles_.material_ids[static_cast<size_t>(particle_index)];
    const PhononMaterial& source = material(old_material);
    const auto old_mode = particles_.modes[static_cast<size_t>(particle_index)];
    const double omega = source.mode_angular_frequency(old_mode);
    const double band_width = 2.0*std::acos(-1.0)*args_.material_interface_frequency_bin_thz;
    const int incident_direction = particles_.velocities[static_cast<size_t>(particle_index)]
        [static_cast<size_t>(args_.material_interface_axis)] > 0.0 ? 1 : 0;

    const double p=args_.material_interface_model==InterfaceModel::AcousticMismatch ? 1.0 :
        args_.material_interface_model==InterfaceModel::MixedMismatch ? args_.material_interface_amm_fraction : 0.0;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const bool use_amm=p>=1 || (p>0 && uniform(thread_rng())<p);
    bool transmit=false;
    int selected_material=old_material,selected_ai=-1;
    if(use_amm) {
        if(std::abs(particles_.velocities[static_cast<size_t>(particle_index)]
            [static_cast<size_t>(args_.material_interface_axis)])<1e-7) {
            particles_.velocities[static_cast<size_t>(particle_index)]
                [static_cast<size_t>(args_.material_interface_axis)]*=-1;
            interface_events_total_.fetch_add(1,std::memory_order_relaxed);
            interface_reflected_total_.fetch_add(1,std::memory_order_relaxed);
            return;
        }
        const auto result=acoustic_interface_sampler_.sample(old_material,next_material,
            source.active_index_for_mode(old_mode),thread_rng());
        selected_material=result.material;selected_ai=result.active;
        transmit=selected_material!=old_material;
    } else {
        // p=0 retains the original DMM draw sequence exactly.
        double source_flux=0,target_flux=0;
        const int reflected_ai=interface_sampler_.sample_band(old_material,1-incident_direction,
            omega,band_width,thread_rng(),source_flux);
        const int transmitted_ai=interface_sampler_.sample_band(next_material,incident_direction,
            omega,band_width,thread_rng(),target_flux);
        const double transmission=target_flux>0 ? target_flux/std::max(1e-300,source_flux+target_flux) : 0;
        transmit=transmitted_ai>=0 && uniform(thread_rng())<transmission;
        selected_material=transmit?next_material:old_material;
        selected_ai=transmit?transmitted_ai:reflected_ai;
    }
    interface_events_total_.fetch_add(1, std::memory_order_relaxed);

    if (selected_ai < 0) {
        if (uses_linearized_transport(args_)) {
            particles_.collision_failed[static_cast<size_t>(particle_index)] = 2;
            collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
            return;
        }
        particles_.velocities[static_cast<size_t>(particle_index)][static_cast<size_t>(args_.material_interface_axis)] *= -1.0;
        interface_reflected_total_.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    const PhononMaterial& destination = material(selected_material);
    const auto new_mode = destination.active_mode_at(selected_ai);
    const double new_omega = destination.mode_angular_frequency(new_mode);
    const double old_deviation = particles_.occupation[static_cast<size_t>(particle_index)] -
        source.bose_occupation(background_temperature_reference(), old_mode);
    const double physical_deviation = material_energy_weight(old_material) * omega * old_deviation;
    const double new_deviation = physical_deviation /
        std::max(1e-300, material_energy_weight(selected_material) * new_omega);
    const double new_occupation = destination.bose_occupation(
        background_temperature_reference(), new_mode) + new_deviation;
    const auto new_velocity = destination.mode_group_velocity(new_mode);
    if (uses_linearized_transport(args_)) {
        const double normal_velocity = new_velocity[static_cast<size_t>(args_.material_interface_axis)];
        const int expected_direction = transmit ? incident_direction : 1 - incident_direction;
        if (!std::isfinite(new_occupation) || !std::isfinite(normal_velocity) ||
            (expected_direction ? !(normal_velocity > 0.0) : !(normal_velocity < 0.0))) {
            particles_.collision_failed[static_cast<size_t>(particle_index)] = 2;
            collision_cache_failures_total_.fetch_add(1, std::memory_order_relaxed);
            return;
        }
    }
    if (!std::isfinite(new_occupation) ||
        (!uses_linearized_transport(args_) && new_occupation < 0.0)) {
        particles_.velocities[static_cast<size_t>(particle_index)][static_cast<size_t>(args_.material_interface_axis)] *= -1.0;
        interface_reflected_total_.fetch_add(1, std::memory_order_relaxed);
        return;
    }
    particles_.material_ids[static_cast<size_t>(particle_index)] = selected_material;
    particles_.modes[static_cast<size_t>(particle_index)] = new_mode;
    particles_.velocities[static_cast<size_t>(particle_index)] = new_velocity;
    particles_.omega[static_cast<size_t>(particle_index)] = new_omega;
    particles_.occupation[static_cast<size_t>(particle_index)] = new_occupation;
    if (transmit) {
        interface_transmitted_total_.fetch_add(1, std::memory_order_relaxed);
    } else {
        interface_reflected_total_.fetch_add(1, std::memory_order_relaxed);
    }
}
