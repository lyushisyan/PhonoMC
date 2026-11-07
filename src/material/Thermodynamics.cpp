#include "PhononMaterial.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// File-free thermodynamic and linewidth queries. Preserve mode reduction order,
// interpolation clamps and monotonic lookup construction, including OpenMP scheduling.
// Modal heat capacity and linearized occupation share the same fixed reference.
double PhononMaterial::mode_heat_capacity(double temperature, const Mode& mode) const {
    const double hw = hbar_ * mode_angular_frequency(mode);
    if (!(temperature > 0.0) || !(hw > 0.0)) {
        return 0.0;
    }
    const double x = hw / (kb_ * temperature);
    if (x < 1e-6) {
        return kb_ * (1.0 - x * x / 12.0);
    }
    if (x > 700.0) {
        return 0.0;
    }
    const double denominator = -std::expm1(-x);
    return kb_ * x * x * std::exp(-x) / (denominator * denominator);
}

double PhononMaterial::linearized_bose_occupation(double temperature,
                                                  double reference_temperature,
                                                  const Mode& mode) const {
    const double hw = hbar_ * mode_angular_frequency(mode);
    if (!(hw > 0.0)) {
        return 0.0;
    }
    return bose_occupation(reference_temperature, mode) +
           mode_heat_capacity(reference_temperature, mode) * (temperature - reference_temperature) /
               hw;
}

// 函数说明：执行线性插值，用于温度与散射参数之间的连续过渡。
double PhononMaterial::lerp(double x0, double x1, double y0, double y1, double x) {
    if (std::abs(x1 - x0) <= 1e-18) {
        return y0;
    }
    const double t = (x - x0) / (x1 - x0);
    return y0 * (1.0 - t) + y1 * t;
}

// 函数说明：在边界钳位条件下做线性插值。
double PhononMaterial::interp_linear_clamped(const std::vector<double>& xs,
                                             const std::vector<double>& ys,
                                             double x) {
    if (xs.empty() || ys.empty()) {
        return 0.0;
    }
    if (xs.size() == 1 || ys.size() == 1) {
        return ys.front();
    }
    if (x <= xs.front()) {
        return ys.front();
    }
    if (x >= xs.back()) {
        return ys.back();
    }
    auto it = std::upper_bound(xs.begin(), xs.end(), x);
    const size_t hi = static_cast<size_t>(std::distance(xs.begin(), it));
    const size_t lo = hi - 1;
    return lerp(xs[lo], xs[hi], ys[lo], ys[hi], x);
}

// 函数说明：计算给定温度下模态的玻色占据数。
double PhononMaterial::bose_occupation(double temperature, const Mode& mode) const {
    const double w = mode_angular_frequency(mode);
    if (temperature <= 0.0 || w <= 0.0) {
        return 0.0;
    }
    const double x = hbar_ * w / (kb_ * temperature);
    if (x > 700.0) {
        return 0.0;
    }
    const double ex = std::exp(x);
    return 1.0 / std::max(ex - 1.0, 1e-12);
}

// 函数说明：计算模态能量（含占据项与零点项）。
double PhononMaterial::mode_energy(double temperature, const Mode& mode) const {
    return hbar_ * mode_angular_frequency(mode) * (bose_occupation(temperature, mode) + 0.5);
}

// 函数说明：按温度和散射率表计算模态寿命。
double PhononMaterial::mode_lifetime(double temperature, const Mode& mode) const {
    const double gamma = mode_gamma(temperature, mode);
    return gamma == 0.0 ? std::numeric_limits<double>::infinity() : (1.0 / (4.0 * M_PI)) / gamma;
}

double
PhononMaterial::mode_relaxation_weight(double temperature, const Mode& mode, double dt) const {
    if (!std::isfinite(dt) || dt < 0.0) {
        throw std::invalid_argument("Relaxation time step must be finite and nonnegative");
    }
    const double gamma = mode_gamma(temperature, mode);
    if (gamma == 0.0 || dt == 0.0) {
        return 0.0;
    }
    return -std::expm1(-(gamma * dt) * (4.0 * M_PI));
}

// Interpolate the finite linewidth, not lifetimes that can be infinite.
double PhononMaterial::mode_gamma(double temperature, const Mode& mode) const {
    if (!std::isfinite(temperature) || temperature < 0.0) {
        throw std::invalid_argument("Lifetime temperature must be finite and nonnegative");
    }
    if (mode[0] < 0 || mode[0] >= qpoint_count_ || mode[1] < 0 || mode[1] >= branch_count_) {
        throw std::out_of_range("Invalid phonon mode for lifetime query");
    }
    const size_t m = static_cast<size_t>(flatten_mode_index(mode));
    const size_t nmode = static_cast<size_t>(qpoint_count_) * static_cast<size_t>(branch_count_);
    if (temperature <= temperature_samples_.front()) {
        return gamma_total_[m];
    }
    if (temperature >= temperature_samples_.back()) {
        return gamma_total_[(temperature_samples_.size() - 1) * nmode + m];
    }
    auto it =
        std::upper_bound(temperature_samples_.begin(), temperature_samples_.end(), temperature);
    const size_t hi = static_cast<size_t>(std::distance(temperature_samples_.begin(), it));
    const size_t lo = hi - 1;
    return lerp(temperature_samples_[lo],
                temperature_samples_[hi],
                gamma_total_[lo * nmode + m],
                gamma_total_[hi * nmode + m],
                temperature);
}

// 函数说明：将模态总量归一化为晶体能量密度量纲。
double PhononMaterial::normalize_to_energy_density(double x) const {
    const double den = std::max(1.0, static_cast<double>(qpoint_count_) * unit_cell_volume_);
    return x / den;
}

// 函数说明：将模态总量归一化为晶体能量密度量纲。
PhononMaterial::Vec3 PhononMaterial::normalize_to_energy_density(const Vec3& x) const {
    const double den = std::max(1.0, static_cast<double>(qpoint_count_) * unit_cell_volume_);
    return {x[0] / den, x[1] / den, x[2] / den};
}

double PhononMaterial::energy_density_normalization() const {
    return std::max(1.0, static_cast<double>(qpoint_count_) * unit_cell_volume_);
}

// 函数说明：累计所有活跃模态得到晶体能量密度。
double PhononMaterial::crystal_energy_density(double temperature) const {
    double e = 0.0;
    for (const Mode& mode : active_mode_list_) {
        const double w = mode_angular_frequency(mode);
        if (w <= 0.0) {
            continue;
        }
        const double n = bose_occupation(temperature, mode);
        e += hbar_ * w * n;
    }
    return normalize_to_energy_density(e) + zero_point_energy_density_;
}

// 函数说明：由温度正向查表得到能量密度（插值），缺表时回退全模态积分。
double PhononMaterial::energy_density_from_temperature(double temperature) const {
    if (temperature_lookup_table_.empty() || energy_lookup_table_.empty()) {
        return crystal_energy_density(temperature);
    }
    return interp_linear_clamped(temperature_lookup_table_, energy_lookup_table_, temperature);
}

// 函数说明：由能量密度反查温度（查找表插值）。
double PhononMaterial::temperature_from_energy_density(double energy_density) const {
    return interp_linear_clamped(energy_lookup_table_, temperature_lookup_table_, energy_density);
}

// 函数说明：构建温度-能量单调查找表用于快速反演。
void PhononMaterial::initialize_temperature_lookup() {
    zero_point_energy_density_ = 0.0;
    for (double w : mode_angular_frequency_data_) {
        zero_point_energy_density_ += 0.5 * hbar_ * std::max(0.0, w);
    }
    zero_point_energy_density_ = normalize_to_energy_density(zero_point_energy_density_);

    // Bose energy depends on the harmonic modes, not the temperatures at which
    // linewidths were computed. In particular, a fixed-rate run at the lowest
    // tabulated temperature must still resolve a slightly colder reservoir.
    const double tmin = 0.0;
    const double tmax =
        temperature_samples_.empty()
            ? 1000.0
            : std::max(1000.0,
                       *std::max_element(temperature_samples_.begin(), temperature_samples_.end()));
    const double dT = std::max(1e-6, temperature_lookup_dt_);
    const double intervals = std::floor((tmax - tmin) / dT + 0.5);
    if (!std::isfinite(intervals) || intervals > std::numeric_limits<int>::max() - 1) {
        throw std::runtime_error("Temperature lookup table exceeds integer indexing capacity");
    }
    const int n = std::max(2, static_cast<int>(intervals) + 1);

    temperature_lookup_table_.assign(static_cast<size_t>(n), 0.0);
    energy_lookup_table_.assign(static_cast<size_t>(n), 0.0);
    std::vector<double> raw_energy(static_cast<size_t>(n), 0.0);

#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < n; ++i) {
        const double T = tmin + static_cast<double>(i) * dT;
        raw_energy[static_cast<size_t>(i)] = crystal_energy_density(T);
    }

    for (int i = 0; i < n; ++i) {
        const double T = tmin + static_cast<double>(i) * dT;
        temperature_lookup_table_[static_cast<size_t>(i)] = T;
        double E = raw_energy[static_cast<size_t>(i)];
        if (i > 0) {
            E = std::max(E,
                         std::nextafter(energy_lookup_table_[static_cast<size_t>(i - 1)],
                                        std::numeric_limits<double>::infinity()));
        }
        energy_lookup_table_[static_cast<size_t>(i)] = E;
    }
}
