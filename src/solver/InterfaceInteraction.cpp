#include "solver/InterfaceInteraction.h"

#include "PhononMaterial.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace phonomc {
void InterfaceInteraction::configure(const InterfaceInteractionConfig& config,
                                     const std::vector<const PhononMaterial*>& materials,
                                     const std::vector<double>& energy_weights) {
    if (!std::isfinite(config.frequency_bin_thz) || config.frequency_bin_thz <= 0) {
        throw std::runtime_error("Interface frequency bin must be finite and positive.");
    }
    if (config.axis < 0 || config.axis > 2 || materials.size() != energy_weights.size()) {
        throw std::invalid_argument("Invalid interface axis or material energy weights.");
    }
    for (size_t i = 0; i < materials.size(); ++i) {
        if (!materials[i] || !std::isfinite(energy_weights[i]) || !(energy_weights[i] > 0)) {
            throw std::invalid_argument("Invalid interface material or energy weight.");
        }
    }
    // Construct a replacement so a rejected configuration cannot mix old/new tables.
    InterfaceInteraction next;
    next.config_ = config;
    next.materials_ = materials;
    next.energy_weights_ = energy_weights;
    next.diffuse_.configure(config.enabled ? materials : std::vector<const PhononMaterial*>{},
                            config.axis);
    if (!std::isfinite(config.amm_fraction) || config.amm_fraction < 0 || config.amm_fraction > 1) {
        throw std::runtime_error("AMM fraction must be in [0,1].");
    }
    if (config.model == InterfaceModel::AcousticMismatch ||
        (config.model == InterfaceModel::MixedMismatch && config.amm_fraction > 0)) {
        next.acoustic_.configure(materials,
                                 config.axis,
                                 config.mass_densities_kg_m3,
                                 config.frequency_bin_thz,
                                 config.parallel_bin_inv_a);
    }
    *this = std::move(next);
}

BoundaryOutcome InterfaceInteraction::interact(const BoundaryState& incoming,
                                               int next_material,
                                               std::mt19937_64& rng) const {
    BoundaryOutcome result;
    result.state = incoming;
    const int old_material = result.state.material;
    const PhononMaterial& source = *materials_.at(static_cast<size_t>(old_material));
    const auto old_mode = result.state.mode;
    const double omega = source.mode_angular_frequency(old_mode);
    const double band_width = 2.0 * std::acos(-1.0) * config_.frequency_bin_thz;
    const int incident_direction =
        result.state.velocity[static_cast<size_t>(config_.axis)] > 0.0 ? 1 : 0;

    double acoustic_probability = 0.0;
    if (config_.model == InterfaceModel::AcousticMismatch) {
        acoustic_probability = 1.0;
    } else if (config_.model == InterfaceModel::MixedMismatch) {
        acoustic_probability = config_.amm_fraction;
    }
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const bool use_amm = acoustic_probability >= 1 ||
                         (acoustic_probability > 0 && uniform(rng) < acoustic_probability);
    bool transmit = false;
    int selected_material = old_material;
    int selected_active_index = -1;
    if (use_amm) {
        if (std::abs(result.state.velocity[static_cast<size_t>(config_.axis)]) < 1e-7) {
            result.state.velocity[static_cast<size_t>(config_.axis)] *= -1;
            ++result.interface_events;
            ++result.interface_reflected;
            result.action = BoundaryAction::Reflect;
            return result;
        }
        const auto selected = acoustic_.sample(
            old_material, next_material, source.active_index_for_mode(old_mode), rng);
        selected_material = selected.material;
        selected_active_index = selected.active;
        transmit = selected_material != old_material;
    } else {
        // Pure DMM skips the model-selection draw. Sample reflection before transmission.
        double source_flux = 0;
        double target_flux = 0;
        const int reflected_active_index = diffuse_.sample_band(
            old_material, 1 - incident_direction, omega, band_width, rng, source_flux);
        const int transmitted_active_index = diffuse_.sample_band(
            next_material, incident_direction, omega, band_width, rng, target_flux);
        const double transmission =
            target_flux > 0 ? target_flux / std::max(1e-300, source_flux + target_flux) : 0;
        transmit = transmitted_active_index >= 0 && uniform(rng) < transmission;
        selected_material = transmit ? next_material : old_material;
        selected_active_index = transmit ? transmitted_active_index : reflected_active_index;
    }
    ++result.interface_events;

    if (selected_active_index < 0) {
        if (config_.linearized) {
            result.error = BoundaryError::InvalidOutgoingMode;
            result.action = BoundaryAction::Fail;
            return result;
        }
        result.state.velocity[static_cast<size_t>(config_.axis)] *= -1.0;
        ++result.interface_reflected;
        result.action = BoundaryAction::Reflect;
        return result;
    }
    const PhononMaterial& destination = *materials_.at(static_cast<size_t>(selected_material));
    const auto new_mode = destination.active_mode_at(selected_active_index);
    const double new_omega = destination.mode_angular_frequency(new_mode);
    const double old_deviation =
        result.state.occupation - source.bose_occupation(config_.reference_temperature, old_mode);
    const double physical_deviation =
        energy_weights_.at(static_cast<size_t>(old_material)) * omega * old_deviation;
    const double new_deviation =
        physical_deviation /
        std::max(1e-300, energy_weights_.at(static_cast<size_t>(selected_material)) * new_omega);
    const double new_occupation =
        destination.bose_occupation(config_.reference_temperature, new_mode) + new_deviation;
    const auto new_velocity = destination.mode_group_velocity(new_mode);
    if (config_.linearized) {
        const double normal_velocity = new_velocity[static_cast<size_t>(config_.axis)];
        const int expected_direction = transmit ? incident_direction : 1 - incident_direction;
        if (!std::isfinite(new_occupation) || !std::isfinite(normal_velocity) ||
            (expected_direction ? !(normal_velocity > 0.0) : !(normal_velocity < 0.0))) {
            result.error = BoundaryError::InvalidOutgoingMode;
            result.action = BoundaryAction::Fail;
            return result;
        }
    }
    if (!std::isfinite(new_occupation) || (!config_.linearized && new_occupation < 0.0)) {
        result.state.velocity[static_cast<size_t>(config_.axis)] *= -1.0;
        ++result.interface_reflected;
        result.action = BoundaryAction::Reflect;
        return result;
    }
    result.state.material = selected_material;
    result.state.mode = new_mode;
    result.state.velocity = new_velocity;
    result.state.omega = new_omega;
    result.state.occupation = new_occupation;
    if (transmit) {
        ++result.interface_transmitted;
        result.action = BoundaryAction::Transmit;
    } else {
        ++result.interface_reflected;
        result.action = BoundaryAction::Reflect;
    }
    return result;
}
} // namespace phonomc
