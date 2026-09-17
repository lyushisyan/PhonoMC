#include "MonteCarloSolver.h"
#include "solver/StageTimer.h"

#include "SimulationDomain.h"
#include <algorithm>
#include <cmath>

void MonteCarloSolver::initialize_local_heat_source(const SimulationDomain& geometry) {
    heat_source_.configure(args_, geometry);
    source_operator_.configure(args_, materials_);
}

double MonteCarloSolver::local_heat_source_integrated_time_factor(
    double time_begin_ps,
    double time_end_ps) const {
    if (!std::isfinite(time_begin_ps) || !std::isfinite(time_end_ps) || time_end_ps <= time_begin_ps) {
        return 0.0;
    }
    const double begin = std::max(time_begin_ps, args_.heat_source_time_start);
    const double end = (args_.heat_source_time_end >= 0.0)
        ? std::min(time_end_ps, args_.heat_source_time_end) : time_end_ps;
    if (end <= begin) {
        return 0.0;
    }
    const double amplitude = std::max(0.0, args_.heat_source_amplitude);
    if (args_.heat_source_time_profile == HeatSourceTimeProfile::Constant) {
        return amplitude * (end - begin);
    }

    const double period = args_.heat_source_period;
    if (!(period > 0.0) || !std::isfinite(period)) {
        return 0.0;
    }
    double on_duration = args_.heat_source_on_duration;
    if (on_duration < 0.0) {
        on_duration = args_.heat_source_duty_cycle * period;
    }
    on_duration = std::clamp(on_duration, 0.0, period);
    if (on_duration <= 0.0) {
        return 0.0;
    }
    const auto accumulated_on_time = [period, on_duration](double relative_time) {
        if (!(relative_time > 0.0)) {
            return 0.0;
        }
        const double cycles = std::floor(relative_time / period);
        const double remainder = relative_time - cycles * period;
        return cycles * on_duration + std::min(remainder, on_duration);
    };
    const double rel_begin = begin - args_.heat_source_time_start;
    const double rel_end = end - args_.heat_source_time_start;
    return amplitude * std::max(0.0, accumulated_on_time(rel_end) - accumulated_on_time(rel_begin));
}

void MonteCarloSolver::apply_local_heat_source_to_occupations(
    double integrated_time_factor_ps) {
    if (!heat_source_.enabled()) return;
    phonomc::StageTimer timer(timer_source_, profile_timers_enabled_);
    heat_source_.accrue(integrated_time_factor_ps);
    const auto deposited = args_.collision_model == CollisionModel::FullMatrix
        ? source_operator_.apply_full_matrix(
            particles_, materials_, {grid_.material_ids, grid_.temperatures, &grid_.cells}, heat_source_,
            particle_spatial_weight_a3_, background_temperature_reference(), rng_)
        : source_operator_.apply(
            particles_, materials_, {grid_.material_ids, grid_.temperatures, &grid_.cells}, heat_source_,
            particle_spatial_weight_a3_, background_temperature_reference());
    if (!deposited.appended_indices.empty()) {
        grid_.cells.invalidate();
        grid_.cells.ensure(particles_, static_cast<int>(grid_.material_ids.size()));
        update_collision_cache(*geometry_, deposited.appended_indices);
    }
    energy_ledger_.record_source(deposited.occupation_updates, deposited.energy_ev);
}
