#include "material/CallawayOperator.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace phonomc {
namespace {
using Vector = std::vector<long double>;
long double dot(const Vector& a, const Vector& b) {
    long double sum = 0;
    for (size_t i=0; i<a.size(); ++i) sum += a[i]*b[i];
    return sum;
}
long double norm(const Vector& a) { return std::sqrt(dot(a,a)); }
std::vector<Vector> gains(const std::vector<double>& rates, std::vector<Vector> basis) {
    std::vector<Vector> orthogonal;
    for (auto v : basis) {
        for (size_t i=0; i<v.size(); ++i) v[i] *= std::sqrt(static_cast<long double>(rates[i]));
        const auto original = norm(v);
        if (original == 0) continue; // inactive mechanism or absent momentum component
        for (auto& x : v) x /= original;
        for (int pass=0; pass<2; ++pass)
            for (const auto& q : orthogonal) {
                const auto projection = dot(q,v);
                for (size_t i=0; i<v.size(); ++i) v[i] -= projection*q[i];
            }
        const auto length = norm(v);
        if (length < 1e-12L) continue; // rank-deficient 1D/2D/support geometry
        for (auto& x : v) x /= length;
        orthogonal.push_back(std::move(v));
    }
    for (auto& q : orthogonal)
        for (size_t i=0; i<q.size(); ++i) q[i] *= std::sqrt(static_cast<long double>(rates[i]));
    return orthogonal;
}
}

CallawayOperator::CallawayOperator(const CallawayData& data) {
    const size_t n=data.omega.size();
    if (!n || data.capacity.size()!=n || data.wavevector.size()!=n ||
        data.normal_rate.size()!=n || data.resistive_rate.size()!=n)
        throw std::invalid_argument("Callaway modal table dimensions differ.");
    scale_.resize(n); rate_.resize(n);
    std::vector<Vector> basis(4,Vector(n));
    long double max_rate=0;
    for (size_t i=0; i<n; ++i) {
        if (!std::isfinite(data.omega[i]) || data.omega[i]<=0 ||
            !std::isfinite(data.capacity[i]) || data.capacity[i]<=0 ||
            !std::isfinite(data.normal_rate[i]) || data.normal_rate[i]<0 ||
            !std::isfinite(data.resistive_rate[i]) || data.resistive_rate[i]<0)
            throw std::invalid_argument("Callaway needs positive frequencies/capacities and finite nonnegative rates.");
        scale_[i]=std::sqrt(static_cast<long double>(data.capacity[i]));
        basis[0][i]=scale_[i];
        for (int a=0; a<3; ++a) {
            if (!std::isfinite(data.wavevector[i][a]))
                throw std::invalid_argument("Nonfinite Callaway wavevector.");
            basis[a+1][i]=scale_[i]*data.wavevector[i][a]/data.omega[i];
        }
        rate_[i]=static_cast<long double>(data.normal_rate[i])+data.resistive_rate[i];
        max_rate=std::max(max_rate,rate_[i]);
    }
    normal_gain_=gains(data.normal_rate,basis);
    resistive_gain_=gains(data.resistive_rate,{basis[0]});
    // Each S is symmetric PSD and bounded above by its rate diagonal.
    norm_bound_=max_rate;
}

std::vector<long double> CallawayOperator::action(const Vector& x) const {
    Vector result(x.size());
    for (size_t i=0; i<x.size(); ++i) result[i]=-rate_[i]*x[i];
    for (const auto* gain : {&normal_gain_, &resistive_gain_})
        for (const auto& u : *gain) {
            const auto coefficient=dot(u,x);
            for (size_t i=0; i<x.size(); ++i) result[i]+=u[i]*coefficient;
        }
    return result;
}

std::vector<double> CallawayOperator::apply(const std::vector<double>& energy) const {
    if (energy.size()!=scale_.size()) throw std::invalid_argument("Callaway energy size mismatch.");
    Vector x(energy.size());
    for (size_t i=0; i<x.size(); ++i) {
        if (!std::isfinite(energy[i])) throw std::invalid_argument("Nonfinite Callaway energy.");
        x[i]=energy[i]/scale_[i];
    }
    const auto y=action(x);
    std::vector<double> out(x.size());
    for (size_t i=0; i<x.size(); ++i) {
        out[i]=static_cast<double>(y[i]*scale_[i]);
        if (!std::isfinite(out[i])) throw std::overflow_error("Callaway action overflow.");
    }
    return out;
}

std::vector<double> CallawayOperator::evolve(const std::vector<double>& energy, double dt) const {
    if (energy.size()!=scale_.size() || !std::isfinite(dt) || dt<0)
        throw std::invalid_argument("Invalid Callaway energy size or time step.");
    Vector x(energy.size());
    long double magnitude=0;
    for (size_t i=0; i<x.size(); ++i) {
        if (!std::isfinite(energy[i])) throw std::invalid_argument("Nonfinite Callaway energy.");
        x[i]=energy[i]/scale_[i];
        magnitude=std::max(magnitude,std::abs(x[i]));
    }
    if (magnitude==0 || dt==0 || norm_bound_==0) return energy;
    // Evaluate the COMPLETE coupled operator; frozen target exponentials do
    // not preserve the collision moments for mode-dependent relaxation times.
    // A scaled norm up to four limits Taylor cancellation to exp(4), while
    // avoiding repeated short expansions for rapidly scattering optical modes.
    const long double steps_ld=std::max(1.0L,std::ceil(dt*norm_bound_/4.0L));
    if (!std::isfinite(steps_ld) || steps_ld>100000 || steps_ld*x.size()>200000000)
        throw std::runtime_error("Callaway exponential work limit exceeded; reduce time_step.");
    const size_t steps=static_cast<size_t>(steps_ld);
    const long double h=static_cast<long double>(dt)/steps;
    const long double scaled_norm=h*norm_bound_;
    for (auto& value : x) value/=magnitude;
    const long double tolerance=2e-14L/steps;
    for (size_t step=0; step<steps; ++step) {
        Vector sum=x, term=x;
        const long double reference=norm(x);
        bool converged=false;
        for (int k=1; k<=40; ++k) {
            term=action(term);
            for (size_t i=0; i<x.size(); ++i) { term[i]*=h/k; sum[i]+=term[i]; }
            // Bound the next term and its remaining geometric majorant using
            // ||hS||_2. The denominator must be positive before testing it.
            if (scaled_norm < k+2) {
                const long double tail=norm(term)*scaled_norm/(k+1)/
                    (1.0L-scaled_norm/(k+2));
                if (tail<=tolerance*reference) { converged=true; break; }
            }
        }
        if (!converged) throw std::runtime_error("Callaway exponential did not converge.");
        x=std::move(sum);
    }
    std::vector<double> result(x.size());
    for (size_t i=0; i<x.size(); ++i) {
        result[i]=static_cast<double>(x[i]*magnitude*scale_[i]);
        if (!std::isfinite(result[i])) throw std::overflow_error("Callaway exponential overflow.");
    }
    return result;
}
} // namespace phonomc
