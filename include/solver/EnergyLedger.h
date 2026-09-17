#pragma once
#include "solver/BoundaryLedger.h"

namespace phonomc {
struct EnergyBalance {
    long long source_updates = 0; // occupation updates, not newly created particles
    double source_energy_ev = 0, lifetime_residual_ev = 0, balance_residual_ev = 0;
    long long drive_updates = 0;
    double drive_energy_ev = 0;
};
// Serial accounting only. All inputs are physical (spatially weighted) eV.
// Signed deviational reservoir energies and signed RTA residuals are intentional.
class EnergyLedger {
public:
    void initialize(double thermal_energy_ev);
    void begin_step();
    void record_source(long long occupation_updates, double energy_ev);
    void record_drive(long long occupation_updates, double energy_ev);
    void set_lifetime_residual(double energy_ev);
    void finish_step(double thermal_energy_ev, const BoundaryBalance& boundary);
    double thermal_energy_ev() const { return thermal_energy_ev_; }
    const EnergyBalance& step() const { return step_; }
    const EnergyBalance& total() const { return total_; }
private:
    bool collecting_ = false;
    double thermal_energy_ev_ = 0;
    EnergyBalance step_, total_;
};
} // namespace phonomc
