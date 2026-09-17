#include "solver/PhononSourceOperator.h"

#include "PhononMaterial.h"
#include "solver/ParticleStorage.h"
#include "solver/PrescribedHeatSource.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace phonomc {
namespace {
constexpr double kHbarEvPs = 6.582119569e-4;
constexpr double kBoltzmannEvK = 8.617333262145e-5;
constexpr double kMinHwEv = 1e-14;
}

void PhononSourceOperator::configure(
    const SimulationConfig& config, const std::vector<const PhononMaterial*>& materials) {
    validate_heat_source_config(config);
    spectrum_ = config.heat_source_spectrum;
    frequency_min_ = config.heat_source_frequency_min;
    frequency_max_ = config.heat_source_frequency_max;
    branches_ = config.heat_source_branches;
    weights_ = config.heat_source_branch_weights;
    linearized_thermal_ = uses_linearized_transport(config);
    gaussian_center_ = config.heat_source_frequency_center.value_or(0);
    gaussian_sigma_ = config.heat_source_frequency_sigma.value_or(1);
    gaussian_fractions_.clear();
    if (!config.heat_source_enabled) return;
    if (!weights_.empty()) {
        const double scale = *std::max_element(weights_.begin(), weights_.end());
        for (double& weight : weights_) weight /= scale;
    }
    for (const auto* material : materials) {
        if (material == nullptr) throw std::invalid_argument("Null source material.");
        for (int branch : branches_) {
            if (branch >= material->branch_count())
                throw std::invalid_argument("heat_source.branches index exceeds a material's branch count.");
        }
        if (spectrum_ == HeatSourceSpectrum::Weighted &&
            weights_.size() != static_cast<size_t>(material->branch_count()))
            throw std::invalid_argument("heat_source.branch_weights must match each material's branch count.");
        if (spectrum_ == HeatSourceSpectrum::Gaussian) {
            const auto& modes = material->active_mode_list();
            size_t positive_modes = 0;
            for (int q = 0; q < material->qpoint_count(); ++q)
                for (int b = 0; b < material->branch_count(); ++b)
                    if (material->mode_angular_frequency({q,b}) > 0) ++positive_modes;
            if (positive_modes != modes.size())
                throw std::invalid_argument("Gaussian source requires all positive-frequency modes, including stationary modes.");
            std::vector<long double> logs(modes.size(), -std::numeric_limits<long double>::infinity());
            long double peak = -std::numeric_limits<long double>::infinity();
            for (size_t m = 0; m < modes.size(); ++m) {
                if (!(mode_weight(*material, modes[m]) > 0)) continue;
                const long double frequency = material->mode_angular_frequency(modes[m]) /
                    6.283185307179586476925286766559L;
                const long double z = (frequency - *config.heat_source_frequency_center) /
                    *config.heat_source_frequency_sigma;
                logs[m] = -0.5L * z * z;
                peak = std::max(peak, logs[m]);
            }
            if (!std::isfinite(peak))
                throw std::invalid_argument("Gaussian source branch/frequency filters select no positive-frequency modes.");
            // Subtract the peak in log space, so a narrow or off-centre Gaussian
            // cannot silently lose its entire power budget to underflow.
            long double sum = 0;
            for (long double& w : logs) { w = std::exp(w - peak); sum += w; }
            std::vector<double> fractions(modes.size());
            for (size_t m = 0; m < modes.size(); ++m) fractions[m] = static_cast<double>(logs[m] / sum);
            gaussian_fractions_.push_back(std::move(fractions));
        }
    }
}

double PhononSourceOperator::mode_weight(
    const PhononMaterial& material, const std::array<int, 2>& mode) const {
    constexpr double two_pi = 6.283185307179586476925286766559;
    const double frequency = material.mode_angular_frequency(mode) / two_pi;
    if (!std::isfinite(frequency) || (frequency_min_ >= 0.0 && frequency < frequency_min_) ||
        (frequency_max_ >= 0.0 && frequency > frequency_max_)) return 0.0;
    if (!branches_.empty() && std::find(branches_.begin(), branches_.end(), mode[1]) == branches_.end())
        return 0.0;
    return spectrum_ == HeatSourceSpectrum::Weighted ? weights_.at(static_cast<size_t>(mode[1])) : 1.0;
}

SourceDeposit PhononSourceOperator::apply(
    ParticleStorage& particles, const std::vector<const PhononMaterial*>& materials,
    const SourceGridView& grid, PrescribedHeatSource& source,
    double particle_volume_a3, double background_temperature) {
    SourceDeposit result;
    if (!source.enabled() || source.pending_energy_ev() <= 0.0) return result;
    assert(particles.aligned());
    const auto& pending = source.cell_pending_ev();
    const size_t ngrid = pending.size();
    if (grid.material_ids.size() != ngrid || grid.temperatures.size() != ngrid)
        throw std::invalid_argument("Source grid and cell-energy ledger sizes disagree.");
    const CellParticleIndex* cells = grid.cells;
    if (cells == nullptr) {
        fallback_cells_.invalidate();
        fallback_cells_.ensure(particles, static_cast<int>(ngrid));
        cells = &fallback_cells_;
    }
    if (!cells->matches(particles, static_cast<int>(ngrid)))
        throw std::invalid_argument("Source requires a current cell index.");
    const auto eligible = [&](int i) {
        const int cell = particles.grid_ids[i];
        if (cell < 0 || cell >= static_cast<int>(ngrid) || pending[cell] <= 0.0 || particles.alive[i] == 0)
            return false;
        if (particles.material_ids[i] != grid.material_ids[cell])
            throw std::runtime_error("Source particle and cell material disagree.");
        const auto& mat = *materials.at(static_cast<size_t>(particles.material_ids[i]));
        const double hw = kHbarEvPs * mat.mode_angular_frequency(particles.modes[i]);
        return std::isfinite(hw) && hw > kMinHwEv && mode_weight(mat, particles.modes[i]) > 0.0;
    };
    for (size_t sv = 0; sv < ngrid; ++sv) {
        if (pending[sv] <= 0.0) continue;
        indices_.clear();
        for (int pos = cells->offsets()[sv]; pos < cells->offsets()[sv + 1]; ++pos) {
            const int i = cells->indices()[pos];
            if (eligible(i)) indices_.push_back(i);
        }
        const int begin = 0, end = static_cast<int>(indices_.size());
        if (end <= begin) continue; // Keep this cell's complete energy budget locally.
        const PhononMaterial& phonon = *materials.at(static_cast<size_t>(grid.material_ids[sv]));
        const double raw_to_physical_ev = static_cast<double>(phonon.active_mode_count()) *
            particles.spatial_volume(grid.material_ids[sv],particle_volume_a3) / phonon.energy_density_normalization();
        if (!(raw_to_physical_ev > 0.0) || !std::isfinite(raw_to_physical_ev))
            throw std::runtime_error("Invalid carrier energy weight for heat source.");
        const double target_raw_energy_ev = pending[sv] / raw_to_physical_ev;
        if (!std::isfinite(target_raw_energy_ev))
            throw std::runtime_error("Non-finite heat-source energy per carrier.");
        increments_.resize(static_cast<size_t>(end - begin));
        if (spectrum_ == HeatSourceSpectrum::Gaussian) {
            // Normalize on the existing sampled quadrature, not per distinct
            // mode. A mode represented twice has twice the quadrature weight.
            std::vector<long double> weights(end);
            long double peak=-std::numeric_limits<long double>::infinity();
            for (int pos=begin;pos<end;++pos) {
                const long double frequency=phonon.mode_angular_frequency(particles.modes[indices_[pos]])/
                    6.283185307179586476925286766559L;
                const long double z=(frequency-gaussian_center_)/gaussian_sigma_;
                weights[pos]=-.5L*z*z;peak=std::max(peak,weights[pos]);
            }
            long double sum=0;
            for(auto& w:weights) {w=std::exp(w-peak);sum+=w;}
            for(int pos=begin;pos<end;++pos) {
                const double hw=kHbarEvPs*phonon.mode_angular_frequency(particles.modes[indices_[pos]]);
                increments_[pos]=static_cast<double>(target_raw_energy_ev*weights[pos]/sum)/hw;
            }
        } else if (spectrum_ == HeatSourceSpectrum::Weighted) {
            double sum = 0.0;
            for (int pos = begin; pos < end; ++pos)
                sum += mode_weight(phonon, particles.modes[indices_[pos]]);
            for (int pos = begin; pos < end; ++pos) {
                const auto& mode = particles.modes[indices_[pos]];
                const double hw = kHbarEvPs * phonon.mode_angular_frequency(mode);
                increments_[pos - begin] = target_raw_energy_ev * (mode_weight(phonon, mode) / sum) / hw;
            }
        } else if (linearized_thermal_) {
            long double capacity=0;
            for(int i:indices_) capacity+=phonon.mode_heat_capacity(background_temperature,particles.modes[i]);
            if (!(capacity>0)) throw std::runtime_error("Source sample has no thermal capacity.");
            for(int pos=begin;pos<end;++pos) {
                const auto mode=particles.modes[indices_[pos]];
                increments_[pos]=static_cast<double>(target_raw_energy_ev*
                    phonon.mode_heat_capacity(background_temperature,mode)/capacity)/
                    (kHbarEvPs*phonon.mode_angular_frequency(mode));
            }
        } else {
            const double base_temperature = std::isfinite(grid.temperatures[sv]) && grid.temperatures[sv] >= 0.0
                ? grid.temperatures[sv] : background_temperature;
            auto thermal_increment_energy = [&](double temperature, double* derivative) {
                double energy = 0.0;
                double slope = 0.0;
                for (int pos = begin; pos < end; ++pos) {
                    const int i = indices_[static_cast<size_t>(pos)];
                    const auto& mode = particles.modes[static_cast<size_t>(i)];
                    const double hw = kHbarEvPs * std::max(0.0, phonon.mode_angular_frequency(mode));
                    const double n0 = phonon.bose_occupation(base_temperature, mode);
                    const double n = phonon.bose_occupation(temperature, mode);
                    energy += hw * (n - n0);
                    if (derivative != nullptr && temperature > 0.0) {
                        const double dndt = hw * n * (n + 1.0) /
                            (kBoltzmannEvK * temperature * temperature);
                        slope += hw * dndt;
                    }
                }
                if (derivative != nullptr) {
                    *derivative = slope;
                }
                return energy;
            };

            double low = std::max(0.0, base_temperature);
            double high = std::max(low + 1.0, 1.05 * low + 1.0);
            double high_energy = thermal_increment_energy(high, nullptr);
            for (int expand = 0; high_energy < target_raw_energy_ev && expand < 60; ++expand) {
                high = 2.0 * high + 1.0;
                high_energy = thermal_increment_energy(high, nullptr);
            }
            if (!(high_energy >= target_raw_energy_ev) || !std::isfinite(high_energy)) {
                throw std::runtime_error("Could not bracket the occupation-source pseudo-temperature.");
            }

            double source_temperature = 0.5 * (low + high);
            const double tolerance = 1e-12 * std::max(1e-6, target_raw_energy_ev);
            for (int iter = 0; iter < 60; ++iter) {
                double derivative = 0.0;
                const double residual =
                    thermal_increment_energy(source_temperature, &derivative) - target_raw_energy_ev;
                if (std::abs(residual) <= tolerance) {
                    break;
                }
                if (residual > 0.0) {
                    high = source_temperature;
                } else {
                    low = source_temperature;
                }
                double candidate = 0.5 * (low + high);
                if (derivative > 0.0 && std::isfinite(derivative)) {
                    const double newton = source_temperature - residual / derivative;
                    if (newton > low && newton < high && std::isfinite(newton)) {
                        candidate = newton;
                    }
                }
                source_temperature = candidate;
            }


            for (int pos = begin; pos < end; ++pos) {
                const auto& mode = particles.modes[indices_[pos]];
                increments_[pos - begin] = phonon.bose_occupation(source_temperature, mode) -
                    phonon.bose_occupation(base_temperature, mode);
            }
        }

        // Stage the cell update, including roundoff correction, before committing.
        double applied = 0.0;
        int correction_pos = begin;
        double largest_energy = -1.0;
        for (int pos = begin; pos < end; ++pos) {
            const double hw = kHbarEvPs * phonon.mode_angular_frequency(particles.modes[indices_[pos]]);
            const double energy = hw * increments_[pos - begin];
            applied += energy;
            if (energy > largest_energy) { largest_energy = energy; correction_pos = pos; }
        }
        const double correction_hw = kHbarEvPs *
            phonon.mode_angular_frequency(particles.modes[indices_[correction_pos]]);
        increments_[correction_pos - begin] += (target_raw_energy_ev - applied) / correction_hw;
        for (int pos = begin; pos < end; ++pos) {
            const double increment = increments_[pos - begin];
            const double next = particles.occupation[indices_[pos]] + increment;
            if (!std::isfinite(increment) || !std::isfinite(next) || next < 0.0)
                throw std::runtime_error("Heat source would produce an invalid occupation; cell budget retained.");
        }
        for (int pos = begin; pos < end; ++pos)
            particles.occupation[indices_[pos]] += increments_[pos - begin];
        result.energy_ev += source.deposit_cell(static_cast<int>(sv));
        result.occupation_updates += end - begin;
    }
    return result;
}

}  // namespace phonomc
