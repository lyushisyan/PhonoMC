#pragma once
#include <array>
#include <vector>

namespace phonomc {
// Fixed-reference, linearized, energy-conserving Callaway model on a supplied
// quadrature: either a full mode grid or the existing equal-weight carriers.
// q is Cartesian angular wavevector (1/Angstrom), omega angular frequency
// (1/ps), capacity eV/K, and rates are inverse lifetimes (1/ps).
// resistive_rate = U + isotope + impurity + defect inverse lifetimes.
struct CallawayData {
    std::vector<double> omega, capacity, normal_rate, resistive_rate;
    std::vector<std::array<double, 3>> wavevector;
};

class CallawayOperator {
public:
    explicit CallawayOperator(const CallawayData& data);
    std::vector<double> apply(const std::vector<double>& energy) const;
    std::vector<double> evolve(const std::vector<double>& energy, double dt_ps) const;
    size_t normal_constraint_rank() const { return normal_gain_.size(); }
private:
    // In x_lambda = energy_lambda / sqrt(capacity_lambda),
    // S_N = R_N - U_N U_N^T, U_N=sqrt(R_N)*orth(sqrt(R_N)*B).
    // B spans energy and crystal momentum; U_R spans only energy.
    std::vector<long double> scale_, rate_;
    std::vector<std::vector<long double>> normal_gain_, resistive_gain_;
    long double norm_bound_ = 0;
    std::vector<long double> action(const std::vector<long double>& x) const;
};
} // namespace phonomc
