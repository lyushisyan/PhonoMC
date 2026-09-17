#include "material/ScatteringMatrix.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace phonomc {
namespace {

[[noreturn]] void invalid(const std::string& reason) {
    throw std::runtime_error("Scattering matrix: " + reason);
}

long double validate_and_norm(const ScatteringMatrixData& data) {
    const std::size_t n = data.modes.size();
    if (!std::isfinite(data.reference_temperature) || data.reference_temperature <= 0.0)
        invalid("reference_temperature must be finite and positive");
    if (n == 0 || n == std::numeric_limits<std::size_t>::max())
        invalid("mode_indices must contain a nonempty, representable mode set");
    if (data.frequency_thz.size() != n) invalid("frequency must have one value per mode");
    std::set<std::array<int, 2>> seen;
    for (std::size_t i = 0; i < n; ++i) {
        if (data.modes[i][0] < 0 || data.modes[i][1] < 0)
            invalid("mode_indices must be nonnegative");
        if (!seen.insert(data.modes[i]).second) invalid("duplicate mode_indices");
        if (!std::isfinite(data.frequency_thz[i]) || data.frequency_thz[i] <= 0.0)
            invalid("frequency must be finite and positive");
    }
    if (data.values.size() != data.row_indices.size())
        invalid("values and row_indices lengths differ");
    if (data.column_offsets.size() != n + 1 || data.column_offsets.front() != 0 ||
        data.column_offsets.back() != data.values.size())
        invalid("column_offsets must span exactly all CSC entries");

    long double norm = 0.0L;
    for (std::size_t col = 0; col < n; ++col) {
        const auto begin = data.column_offsets[col], end = data.column_offsets[col + 1];
        if (end < begin || end > data.values.size()) invalid("column_offsets must be nondecreasing and in bounds");
        long double sum = 0.0L, abs_sum = 0.0L, diagonal = 0.0L;
        for (std::size_t p = begin; p < end; ++p) {
            const auto row = data.row_indices[p];
            if (row >= n) invalid("row_indices entry outside mode range");
            if (p > begin && row <= data.row_indices[p - 1])
                invalid("row_indices must be strictly increasing within each column (no duplicates)");
            const double value = data.values[p];
            if (!std::isfinite(value)) invalid("values must be finite");
            sum += static_cast<long double>(value);
            abs_sum += std::abs(static_cast<long double>(value));
            if (row == col) diagonal = value;
        }
        if (!std::isfinite(abs_sum)) invalid("column norm exceeds numeric range");
        if (std::abs(sum) > 1.0e-9L * abs_sum)
            invalid("column " + std::to_string(col) + " does not conserve energy (column sum is not zero)");
        if (diagonal > 1.0e-12L * abs_sum)
            invalid("diagonal must be nonpositive");
        norm = std::max(norm, abs_sum);
    }
    return norm;
}

long double norm_one(const std::vector<long double>& vector) {
    long double norm = 0.0L;
    for (const auto value : vector) norm += std::abs(value);
    if (!std::isfinite(norm)) invalid("exponential action produced a non-finite intermediate result");
    return norm;
}

} // namespace

void validate_scattering_matrix_data(const ScatteringMatrixData& data) {
    (void)validate_and_norm(data);
}

ScatteringMatrix::ScatteringMatrix(ScatteringMatrixData data) : data_(std::move(data)) {
    norm_one_ = validate_and_norm(data_);
}

std::vector<double> ScatteringMatrix::evolve(std::vector<double> energy, double dt) const {
    if (energy.size() != size()) invalid("energy vector size does not match modes");
    if (!std::isfinite(dt) || dt < 0.0) invalid("dt must be finite and nonnegative");
    double energy_scale = 0.0;
    for (double value : energy) {
        if (!std::isfinite(value)) invalid("energy vector must be finite");
        energy_scale = std::max(energy_scale, std::abs(value));
    }
    if (dt == 0.0 || norm_one_ == 0.0L || energy_scale == 0.0) return energy;

    // Each substep has ||h A||_1 <= 1, so its Taylor tail is bounded by a
    // decreasing geometric series. Scaling the input avoids large-energy
    // intermediate overflow; accumulation in long double limits stiff-step
    // roundoff without changing the operator or projecting conserved modes.
    const long double scaled_norm = norm_one_ * static_cast<long double>(dt);
    constexpr std::size_t max_steps = 100000;
    constexpr std::size_t max_terms = 64;
    constexpr std::size_t max_sparse_work = 2000000000;
    if (!std::isfinite(scaled_norm) || scaled_norm > max_steps)
        invalid("exponential action work limit exceeded; reduce dt or matrix rates");
    const auto steps = static_cast<std::size_t>(std::max(1.0L, std::ceil(scaled_norm)));
    // Conservative bound includes vector work for matrices with empty columns.
    const long double work = static_cast<long double>(steps) * max_terms *
        (static_cast<long double>(data_.values.size()) + static_cast<long double>(size()));
    if (work > max_sparse_work)
        invalid("exponential action work limit exceeded; reduce dt or matrix size");
    const long double h = static_cast<long double>(dt) / steps;
    const long double h_norm = scaled_norm / steps;
    const long double relative_tolerance = std::max(
        8.0L * std::numeric_limits<long double>::epsilon(),
        0.125L * std::numeric_limits<double>::epsilon() / steps);

    std::vector<long double> state(size()), term(size()), sum(size()), product(size());
    for (std::size_t i = 0; i < size(); ++i)
        state[i] = static_cast<long double>(energy[i]) / energy_scale;
    for (std::size_t step = 0; step < steps; ++step) {
        term = state;
        sum = state;
        bool converged = false;
        for (std::size_t k = 1; k <= max_terms; ++k) {
            std::fill(product.begin(), product.end(), 0.0L);
            const long double factor = h / k;
            for (std::size_t col = 0; col < size(); ++col) {
                if (term[col] == 0.0L) continue;
                for (std::size_t p = data_.column_offsets[col]; p < data_.column_offsets[col + 1]; ++p)
                    product[data_.row_indices[p]] +=
                        (factor * static_cast<long double>(data_.values[p])) * term[col];
            }
            term.swap(product);
            for (std::size_t i = 0; i < size(); ++i) sum[i] += term[i];
            const long double term_norm = norm_one(term);
            const long double sum_norm = norm_one(sum);
            const long double remainder_bound = term_norm * h_norm / (k + 1) /
                (1.0L - h_norm / (k + 2));
            if (term_norm == 0.0L || remainder_bound <= relative_tolerance * sum_norm) {
                converged = true;
                break;
            }
        }
        if (!converged) invalid("exponential action Taylor series did not converge");
        state.swap(sum);
    }
    for (std::size_t i = 0; i < size(); ++i) {
        const long double value = state[i] * energy_scale;
        if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<double>::max())
            invalid("exponential action result exceeds finite energy range");
        energy[i] = static_cast<double>(value);
    }
    return energy;
}

} // namespace phonomc
