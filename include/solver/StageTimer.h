#pragma once
#include <chrono>

namespace phonomc {
// Exclusive stage scope; nested calls must not charge the same elapsed time.
class StageTimer {
public:
    StageTimer(double& accumulator, bool enabled) : accumulator_(accumulator), enabled_(enabled) {
        if (enabled_) begin_ = Clock::now();
    }
    ~StageTimer() {
        if (enabled_) accumulator_ += std::chrono::duration<double>(Clock::now() - begin_).count();
    }
    StageTimer(const StageTimer&) = delete;
    StageTimer& operator=(const StageTimer&) = delete;
private:
    using Clock = std::chrono::steady_clock;
    double& accumulator_;
    bool enabled_;
    Clock::time_point begin_;
};
}  // namespace phonomc
