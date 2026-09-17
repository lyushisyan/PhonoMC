#pragma once
#include <array>
#include <cstddef>
#include <vector>

namespace phonomc {

// Linearized mode ENERGY generator: dE/dt = A E, with time in ps.
// Columns and rows share the explicit (zero-based qpoint, branch) mode order.
// CSC entries are sorted by row within each column. Signed off-diagonals are
// allowed: this is a linearized collision operator, not a Markov generator.
struct ScatteringMatrixData {
    double reference_temperature = 0.0;
    std::vector<std::array<int, 2>> modes;
    std::vector<double> frequency_thz;
    std::vector<std::size_t> column_offsets;
    std::vector<std::size_t> row_indices;
    std::vector<double> values;
};

// Checks structure, finite data, diagonal signs and energy conservation.
// Thermal equilibrium and agreement with a material are checked when bound.
void validate_scattering_matrix_data(const ScatteringMatrixData& data);

class ScatteringMatrix {
public:
    explicit ScatteringMatrix(ScatteringMatrixData data);
    const ScatteringMatrixData& data() const noexcept { return data_; }
    std::size_t size() const noexcept { return data_.modes.size(); }

    // Sparse exponential action exp(dt A) E. No dense propagator is formed,
    // and signed energies are never clipped. Throws on invalid input,
    // non-finite results, lack of convergence, or excessive requested work.
    std::vector<double> evolve(std::vector<double> energy, double dt) const;

private:
    ScatteringMatrixData data_;
    long double norm_one_ = 0.0L;
};

} // namespace phonomc
