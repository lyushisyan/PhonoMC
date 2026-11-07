#pragma once
#include "solver/RoughBoundarySampler.h"

namespace phonomc {
// Value input/output for one event. No storage references survive compaction.
struct BoundaryState {
    std::array<int, 2> mode{};
    int material = 0;
    std::array<double, 3> velocity{};
    double omega = 0;
    double occupation = 0;
};
enum class BoundaryAction { Unchanged, Absorb, Periodic, Reflect, Transmit, Fail };
enum class BoundaryError { None, InvalidOutgoingMode };

struct BoundaryOutcome {
    BoundaryAction action = BoundaryAction::Unchanged;
    BoundaryState state;
    std::array<double, 3> position_shift{};
    double absorbed_energy_ev = 0;
    bool clear_cached_energy = false;
    BoundaryError error = BoundaryError::None;
    RoughSampleEvents rough;
    long long interface_events = 0;
    long long interface_transmitted = 0;
    long long interface_reflected = 0;
};
} // namespace phonomc
