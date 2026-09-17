#pragma once
#include "SimulationConfig.h"
#include <algorithm>
#include <cmath>
#include <deque>
#include <stdexcept>
#include <vector>

namespace phonomc {
// Non-overlapping time windows of global and layer conductivities. Passing this
// stationarity screen is not an estimate of error relative to the infinite-time limit.
class ConvergenceMonitor {
public:
    explicit ConvergenceMonitor(const SimulationConfig& c = {}) : config_(c), end_(c.convergence_min_time_ps+c.convergence_window_ps) {}
    void sample(double t,const std::vector<double>& values) {
        if (!config_.convergence_stop || converged_) return;
        if (!std::isfinite(t) || values.empty() ||
            std::any_of(values.begin(),values.end(),[](double v){return !std::isfinite(v);}))
            throw std::invalid_argument("Non-finite convergence sample.");
        if (previous_.empty()) {previous_=values;previous_time_=t;integral_.assign(values.size(),0);return;}
        if (t<=previous_time_ || values.size()!=previous_.size()) throw std::invalid_argument("Invalid convergence sample order/shape.");
        double a=std::max(previous_time_,config_.convergence_min_time_ps);
        while (a<t && !converged_) {
            const double b=std::min(t,end_);
            for (size_t i=0;i<values.size();++i) {
                const double slope=(values[i]-previous_[i])/(t-previous_time_);
                integral_[i]+=(b-a)*(previous_[i]+slope*((a+b)/2-previous_time_));
            }
            a=b;
            if (b>=end_-1e-10) {
                std::vector<double> mean;
                for (auto v:integral_) mean.push_back(static_cast<double>(v/config_.convergence_window_ps));
                blocks_.push_back(mean);++completed_windows_;
                if (blocks_.size()>static_cast<size_t>(config_.convergence_windows)) blocks_.pop_front();
                last_end_=end_;end_+=config_.convergence_window_ps;
                std::fill(integral_.begin(),integral_.end(),0);
                evaluate();
            }
        }
        previous_=values;previous_time_=t;
    }
    bool converged() const {return converged_;}
    int completed_windows() const {return completed_windows_;}
    double last_window_end() const {return last_end_;}
    double maximum_spread() const {return spread_;}
    const std::deque<std::vector<double>>& blocks() const {return blocks_;}
private:
    void evaluate() {
        if (blocks_.size()<static_cast<size_t>(config_.convergence_windows)) return;
        double average=0;for (const auto& b:blocks_) average+=b[0]/blocks_.size();
        // Zero current / a ballistic initial plateau must not qualify as transport convergence.
        if (average<=config_.convergence_absolute_tolerance) return;
        spread_=0;
        for (size_t i=0;i<blocks_.front().size();++i) {
            double lo=blocks_.front()[i],hi=lo;
            for (const auto& b:blocks_) {lo=std::min(lo,b[i]);hi=std::max(hi,b[i]);}
            spread_=std::max(spread_,hi-lo);
        }
        converged_=spread_<=config_.convergence_absolute_tolerance+config_.convergence_relative_tolerance*std::abs(average);
    }
    SimulationConfig config_;
    bool converged_=false;
    int completed_windows_=0;
    double previous_time_=0,end_=0,last_end_=0,spread_=0;
    std::vector<double> previous_;
    std::vector<long double> integral_;
    std::deque<std::vector<double>> blocks_;
};
}
