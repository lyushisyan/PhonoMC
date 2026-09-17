#pragma once

#include <vector>
#include <array>

struct SimulationConfig;
class SimulationDomain;

namespace phonomc {

// Prescribed lattice power, not an electron transport or electron-phonon model.
// Owns a separate energy reservoir for every cell. Energy enters the thermal
// state only when deposited. The linearized source can seed empty cells;
// the legacy RTA source defers their budget until carriers are available.
class PrescribedHeatSource {
public:
    void configure(const SimulationConfig& config, const SimulationDomain& geometry);
    void accrue(double integrated_time_ps);
    double deposit_cell(int cell);
    bool enabled() const { return enabled_; }
    double base_power_w() const;
    double prescribed_energy_ev() const;
    double deposited_energy_ev() const;
    double pending_energy_ev() const;
    const std::vector<double>& cell_rates_evps() const { return rates_; }
    const std::vector<double>& cell_prescribed_ev() const { return prescribed_; }
    const std::vector<double>& cell_deposited_ev() const { return deposited_; }
    const std::vector<double>& cell_pending_ev() const { return pending_; }
    const std::array<double,3>& cell_source_position(size_t cell) const { return positions_.at(cell); }

private:
    bool enabled_ = false;
    std::vector<double> rates_, prescribed_, deposited_, pending_;
    std::vector<std::array<double,3>> positions_;
};

}  // namespace phonomc
