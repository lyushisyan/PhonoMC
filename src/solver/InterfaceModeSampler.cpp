#include "solver/InterfaceModeSampler.h"
#include "PhononMaterial.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace phonomc {
void InterfaceModeSampler::configure(const std::vector<const PhononMaterial*>& materials, int axis) {
    if (axis < 0 || axis > 2) throw std::invalid_argument("Interface sampling axis must be 0, 1 or 2.");
    std::vector<std::array<Bank, 2>> banks(materials.size());
    for (int material_id = 0; material_id < static_cast<int>(materials.size()); ++material_id) {
        if (materials[material_id] == nullptr) throw std::invalid_argument("Null interface material.");
        const PhononMaterial& mat = *materials[material_id];
        for (int direction = 0; direction < 2; ++direction) {
            auto& bank = banks[static_cast<size_t>(material_id)][static_cast<size_t>(direction)];
            std::vector<std::pair<double, int>> ordered;
            for (int ai = 0; ai < mat.active_mode_count(); ++ai) {
                const auto mode = mat.active_mode_at(ai);
                const double component = mat.mode_group_velocity(mode)[static_cast<size_t>(axis)];
                if ((direction == 0 && component < 0.0) || (direction == 1 && component > 0.0)) {
                    ordered.push_back({mat.mode_angular_frequency(mode), ai});
                }
            }
            std::sort(ordered.begin(), ordered.end());
            bank.flux_prefix.assign(ordered.size() + 1, 0.0);
            for (size_t i = 0; i < ordered.size(); ++i) {
                bank.angular_frequencies.push_back(ordered[i].first);
                bank.active_indices.push_back(ordered[i].second);
                const auto mode = mat.active_mode_at(ordered[i].second);
                const double flux = std::abs(mat.mode_group_velocity(mode)[static_cast<size_t>(axis)]) /
                    mat.energy_density_normalization();
                bank.flux_prefix[i + 1] = bank.flux_prefix[i] + flux;
            }
        }
    }
    banks_.swap(banks);
}

int InterfaceModeSampler::sample(
    int material_id,
    int direction,
    double omega,
    double tolerance,
    std::mt19937_64& rng,
    double& flux_sum) const {
    if (!std::isfinite(omega) || !std::isfinite(tolerance) || tolerance < 0)
        throw std::invalid_argument("Interface frequency must be finite; tolerance must be finite and nonnegative.");
    const auto& bank = banks_.at(static_cast<size_t>(material_id)).at(static_cast<size_t>(direction));
    const auto first = std::lower_bound(bank.angular_frequencies.begin(), bank.angular_frequencies.end(), omega - tolerance);
    const auto last = std::upper_bound(bank.angular_frequencies.begin(), bank.angular_frequencies.end(), omega + tolerance);
    const size_t begin = static_cast<size_t>(first - bank.angular_frequencies.begin());
    const size_t end = static_cast<size_t>(last - bank.angular_frequencies.begin());
    return draw(bank,begin,end,rng,flux_sum);
}

int InterfaceModeSampler::sample_band(int material_id, int direction, double omega,
    double band_width, std::mt19937_64& rng, double& flux_sum) const {
    if (!std::isfinite(omega) || omega<0 || !std::isfinite(band_width) || band_width<=0)
        throw std::invalid_argument("Invalid shared interface frequency band.");
    const auto& bank=banks_.at(static_cast<size_t>(material_id)).at(static_cast<size_t>(direction));
    const double bin=std::floor(omega/band_width);
    const double lo=bin*band_width, hi=(bin+1)*band_width;
    const auto begin=std::lower_bound(bank.angular_frequencies.begin(),bank.angular_frequencies.end(),lo);
    const auto end=std::lower_bound(bank.angular_frequencies.begin(),bank.angular_frequencies.end(),hi);
    return draw(bank,begin-bank.angular_frequencies.begin(),end-bank.angular_frequencies.begin(),rng,flux_sum);
}

int InterfaceModeSampler::draw(const Bank& bank,size_t begin,size_t end,
    std::mt19937_64& rng,double& flux_sum) const {
    flux_sum = (end > begin) ? bank.flux_prefix[end] - bank.flux_prefix[begin] : 0.0;
    if (!(flux_sum > 0.0)) {
        return -1;
    }
    std::uniform_real_distribution<double> draw(0.0, flux_sum);
    const double target = bank.flux_prefix[begin] + draw(rng);
    const auto selected = std::upper_bound(
        bank.flux_prefix.begin() + static_cast<std::ptrdiff_t>(begin + 1),
        bank.flux_prefix.begin() + static_cast<std::ptrdiff_t>(end + 1),
        target);
    size_t position = static_cast<size_t>(selected - bank.flux_prefix.begin());
    position = std::clamp(position, begin + 1, end) - 1;
    return bank.active_indices[position];
}

}  // namespace phonomc
