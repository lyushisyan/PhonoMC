#include "solver/BoundaryLedger.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace phonomc {
void BoundaryLedger::configure(const std::vector<int>& counts) {
    if (phase_ == Phase::Collecting || phase_ == Phase::Reduced)
        throw std::logic_error("Cannot reconfigure an unfinished boundary step.");
    if (std::any_of(counts.begin(), counts.end(), [](int n) { return n < 0; }))
        throw std::invalid_argument("Negative initial reservoir emission count.");
    auto previous = counts;
    std::vector<int> current(counts.size(), 0);
    previous_.swap(previous);
    current_.swap(current);
    workers_.clear();
    leaving_by_worker_.clear();
    step_ = {};
    total_ = {};
    invalid_event_.store(false, std::memory_order_relaxed);
    phase_ = Phase::Ready;
}
void BoundaryLedger::begin_step(int workers) {
    if (phase_ != Phase::Ready && phase_ != Phase::Committed)
        throw std::logic_error("Previous boundary step has not been committed.");
    if (workers < 1) throw std::invalid_argument("Boundary worker count must be positive.");
    const auto count = static_cast<size_t>(workers);
    if (!previous_.empty() && count > leaving_by_worker_.max_size() / previous_.size())
        throw std::length_error("Boundary worker/reservoir dimensions overflow.");
    // Reuse allocations when the worker count and reservoir layout are unchanged.
    workers_.resize(count);
    leaving_by_worker_.resize(count * previous_.size());
    std::fill(workers_.begin(), workers_.end(), Worker{});
    std::fill(leaving_by_worker_.begin(), leaving_by_worker_.end(), 0);
    std::fill(current_.begin(), current_.end(), 0);
    step_ = {};
    invalid_event_.store(false, std::memory_order_relaxed);
    phase_ = Phase::Collecting;
}
void BoundaryLedger::record_absorption(int worker, int reservoir, double energy_ev) noexcept {
    if (phase_ != Phase::Collecting || worker < 0 || static_cast<size_t>(worker) >= workers_.size() ||
        reservoir < -1 || (reservoir >= 0 && static_cast<size_t>(reservoir) >= previous_.size()) ||
        !std::isfinite(energy_ev)) {
        invalid_event_.store(true, std::memory_order_relaxed);
        return;
    }
    // -1 means an absorbing facet with no refill reservoir, not an invalid event.
    if (reservoir >= 0) {
        auto& leaving = leaving_by_worker_[static_cast<size_t>(worker) * previous_.size() + reservoir];
        if (leaving == std::numeric_limits<int>::max()) {
            invalid_event_.store(true, std::memory_order_relaxed);
            return;
        }
        ++leaving;
    }
    auto& slot = workers_[worker];
    ++slot.absorbed;
    slot.energy_ev += energy_ev;
}
void BoundaryLedger::record_injection(double energy_ev) {
    if (phase_ != Phase::Collecting) throw std::logic_error("Injection outside boundary collection.");
    if (!std::isfinite(energy_ev)) throw std::invalid_argument("Non-finite injected energy.");
    ++step_.injected;
    step_.injected_energy_ev += energy_ev;
}
void BoundaryLedger::finish_transport() {
    if (phase_ != Phase::Collecting) throw std::logic_error("Boundary transport already reduced or not started.");
    if (invalid_event_.load(std::memory_order_relaxed))
        throw std::runtime_error("Invalid boundary event (worker, reservoir, energy or count overflow).");
    BoundaryBalance reduced = step_;
    // Ascending worker order matches the old OpenMP energy reduction exactly.
    for (const auto& slot : workers_) {
        reduced.absorbed += slot.absorbed;
        reduced.absorbed_energy_ev += slot.energy_ev;
    }
    if (!std::isfinite(reduced.absorbed_energy_ev) || !std::isfinite(reduced.injected_energy_ev))
        throw std::overflow_error("Boundary energy reduction overflow.");
    for (size_t r = 0; r < previous_.size(); ++r) {
        long long leaving = 0;
        for (size_t worker = 0; worker < workers_.size(); ++worker)
            leaving += leaving_by_worker_[worker * previous_.size() + r];
        if (leaving > std::numeric_limits<int>::max())
            throw std::overflow_error("Reservoir refill count exceeds int capacity.");
        current_[r] = static_cast<int>(leaving);
    }
    step_ = reduced;
    phase_ = Phase::Reduced;
}
void BoundaryLedger::commit_step() {
    if (phase_ != Phase::Reduced) throw std::logic_error("Boundary commit requires one completed reduction.");
    if (invalid_event_.load(std::memory_order_relaxed))
        throw std::runtime_error("Boundary event recorded after transport reduction.");
    if (step_.absorbed > std::numeric_limits<long long>::max() - total_.absorbed ||
        step_.injected > std::numeric_limits<long long>::max() - total_.injected ||
        !std::isfinite(total_.absorbed_energy_ev + step_.absorbed_energy_ev) ||
        !std::isfinite(total_.injected_energy_ev + step_.injected_energy_ev))
        throw std::overflow_error("Cumulative boundary balance overflow.");
    total_.absorbed += step_.absorbed;
    total_.injected += step_.injected;
    total_.absorbed_energy_ev += step_.absorbed_energy_ev;
    total_.injected_energy_ev += step_.injected_energy_ev;
    previous_ = current_;
    phase_ = Phase::Committed;
}
int BoundaryLedger::emission_count(int reservoir) const {
    return previous_.at(static_cast<size_t>(reservoir));
}
} // namespace phonomc
