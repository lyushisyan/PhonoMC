#pragma once
#include <atomic>
#include <vector>

namespace phonomc {
struct BoundaryBalance {
    long long absorbed = 0, injected = 0;
    double absorbed_energy_ev = 0, injected_energy_ev = 0;
    long long net() const { return injected - absorbed; }
};

// No particles, geometry or OpenMP dependency. Each worker slot has one writer.
// Except record_absorption(), calls run only outside parallel regions. A failed
// step is not resumable; configure() is not a recovery mechanism for that step.
class BoundaryLedger {
public:
    void configure(const std::vector<int>& initial_emission_counts);
    void begin_step(int workers);
    void record_absorption(int worker, int reservoir, double energy_ev) noexcept;
    void record_injection(double energy_ev);
    void finish_transport();
    void commit_step();
    int emission_count(int reservoir) const;
    const BoundaryBalance& step() const { return step_; }
    const BoundaryBalance& total() const { return total_; }
private:
    enum class Phase { Ready, Collecting, Reduced, Committed };
    Phase phase_ = Phase::Ready;
    struct Worker {
        long long absorbed = 0;
        double energy_ev = 0;
    };
    std::vector<Worker> workers_;
    std::vector<int> leaving_by_worker_, previous_, current_;
    std::atomic<bool> invalid_event_ {false};
    BoundaryBalance step_, total_;
};
} // namespace phonomc
