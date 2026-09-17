#include "solver/EnergyLedger.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace phonomc {
namespace {
void finite(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("Non-finite energy ledger input.");
}
}
void EnergyLedger::initialize(double energy) {
    if (collecting_) throw std::logic_error("Cannot reset an unfinished energy step.");
    finite(energy);
    thermal_energy_ev_ = energy;
    step_ = {};
    total_ = {};
}
void EnergyLedger::begin_step() {
    if (collecting_) throw std::logic_error("Previous energy step has not finished.");
    step_ = {};
    collecting_ = true;
}
void EnergyLedger::record_source(long long updates, double energy) {
    if (!collecting_) throw std::logic_error("Source accounting outside energy step.");
    if (updates < 0) throw std::invalid_argument("Negative source occupation-update count.");
    finite(energy);
    if (updates > std::numeric_limits<long long>::max() - step_.source_updates)
        throw std::overflow_error("Source occupation-update count overflow.");
    finite(step_.source_energy_ev + energy);
    step_.source_updates += updates;
    step_.source_energy_ev += energy;
}
void EnergyLedger::set_lifetime_residual(double energy) {
    if (!collecting_) throw std::logic_error("RTA accounting outside energy step.");
    finite(energy);
    step_.lifetime_residual_ev = energy;
}
void EnergyLedger::finish_step(double energy, const BoundaryBalance& boundary) {
    if (!collecting_) throw std::logic_error("Energy step already finished or not started.");
    finite(energy);
    finite(boundary.injected_energy_ev);
    finite(boundary.absorbed_energy_ev);
    // Preserve the expression's evaluation order: this is accounting, not a
    // correction to occupations or to any of the supplied exchange energies.
    const double residual = energy - thermal_energy_ev_ - step_.source_energy_ev - step_.drive_energy_ev -
        boundary.injected_energy_ev + boundary.absorbed_energy_ev - step_.lifetime_residual_ev;
    finite(residual);
    if (step_.source_updates > std::numeric_limits<long long>::max() - total_.source_updates)
        throw std::overflow_error("Cumulative source occupation-update count overflow.");
    finite(total_.source_energy_ev + step_.source_energy_ev);
    finite(total_.lifetime_residual_ev + step_.lifetime_residual_ev);
    finite(total_.balance_residual_ev + residual);
    finite(total_.drive_energy_ev + step_.drive_energy_ev);
    if (step_.drive_updates>std::numeric_limits<long long>::max()-total_.drive_updates)
        throw std::overflow_error("Cumulative gradient drive update count overflow.");
    step_.balance_residual_ev = residual;
    total_.source_updates += step_.source_updates;
    total_.source_energy_ev += step_.source_energy_ev;
    total_.lifetime_residual_ev += step_.lifetime_residual_ev;
    total_.balance_residual_ev += step_.balance_residual_ev;
    total_.drive_updates += step_.drive_updates;
    total_.drive_energy_ev += step_.drive_energy_ev;
    thermal_energy_ev_ = energy;
    collecting_ = false;
}
void EnergyLedger::record_drive(long long updates, double energy) {
    if (!collecting_) throw std::logic_error("Gradient drive accounting outside energy step.");
    if (updates<0 || updates>std::numeric_limits<long long>::max()-step_.drive_updates)
        throw std::overflow_error("Invalid gradient drive update count.");
    finite(energy);finite(step_.drive_energy_ev+energy);
    step_.drive_updates+=updates;step_.drive_energy_ev+=energy;
}
} // namespace phonomc
