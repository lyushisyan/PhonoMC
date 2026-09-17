#include "PhononMaterial.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

phonomc::CallawayData PhononMaterial::linearized_rta_data(double temperature) const {
    if (!std::isfinite(temperature) || temperature <= 0 || temperature_samples_.empty() ||
        temperature < temperature_samples_.front() || temperature > temperature_samples_.back())
        throw std::invalid_argument("Linearized RTA reference temperature is outside the scattering data range.");
    phonomc::CallawayData data;
    for (const auto& mode : active_mode_list_) {
        data.omega.push_back(mode_angular_frequency(mode));
        data.capacity.push_back(mode_heat_capacity(temperature, mode));
        data.wavevector.push_back({0,0,0});
        data.normal_rate.push_back(0);
        data.resistive_rate.push_back(4*std::acos(-1.0)*mode_gamma(temperature, mode));
    }
    return data;
}

PhononMaterial::Vec3 PhononMaterial::mode_wavevector(const Mode& mode) const {
    if (mode[0]<0 || mode[0]>=qpoint_count_ || mode[1]<0 || mode[1]>=branch_count_)
        throw std::out_of_range("Invalid wavevector mode.");
    return mode_wavevector_data_.at(static_cast<size_t>(flatten_mode_index(mode)));
}

phonomc::CallawayData PhononMaterial::callaway_data(double temperature) const {
    if (gamma_normal_.empty() || gamma_umklapp_.empty())
        throw std::invalid_argument("Callaway requires gamma_N and gamma_U in kappa-fbz.hdf5; regenerate with phono3py --br --nu and expand consistently. No N/U split is inferred.");
    if (!std::isfinite(temperature) || temperature<=0 || temperature<temperature_samples_.front() ||
        temperature>temperature_samples_.back())
        throw std::invalid_argument("Callaway reference temperature is outside the N/U data range.");
    for (size_t i=0; i<gamma_table_.size(); ++i) {
        const double sum=gamma_normal_[i]+gamma_umklapp_[i];
        if (!std::isfinite(sum) || std::abs(sum-gamma_table_[i])>1e-7*std::max(sum,gamma_table_[i])+1e-14)
            throw std::invalid_argument("Callaway data inconsistent: gamma_N + gamma_U != gamma. Supply unmodified, same-run linewidths in the same units.");
    }
    // Check the actual Gamma-centred mesh, not just its point count. Retain
    // original first-BZ representatives: independent component wrapping is
    // incorrect for a nonorthogonal reciprocal lattice.
    std::set<std::array<int,3>> grid_points;
    for (const auto& q : qpoint_fractions_) {
        std::array<int,3> key;
        for (int a=0; a<3; ++a) {
            const double fractional=q[a]-std::floor(q[a]);
            const double scaled=fractional*mesh_[a];
            const long index=std::lround(scaled);
            if (std::abs(scaled-index)>2e-6*mesh_[a])
                throw std::invalid_argument("Callaway requires a complete Gamma-centred full q mesh.");
            key[a]=static_cast<int>(index%mesh_[a]);
        }
        if (!grid_points.insert(key).second)
            throw std::invalid_argument("Callaway q mesh contains duplicate points modulo reciprocal translations.");
    }
    // Bound all reciprocal translations that could shorten each q. The inverse
    // row norm bounds each fractional component of a Cartesian norm ball.
    const auto& b=reciprocal_;
    const double det=b[0][0]*(b[1][1]*b[2][2]-b[1][2]*b[2][1])-
        b[0][1]*(b[1][0]*b[2][2]-b[1][2]*b[2][0])+b[0][2]*(b[1][0]*b[2][1]-b[1][1]*b[2][0]);
    if (!std::isfinite(det) || det==0) throw std::invalid_argument("Callaway reciprocal lattice is singular.");
    double inv[3][3];
    for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
        inv[i][j]=(b[(j+1)%3][(i+1)%3]*b[(j+2)%3][(i+2)%3]-
                   b[(j+1)%3][(i+2)%3]*b[(j+2)%3][(i+1)%3])/det;
    for (int qi=0; qi<qpoint_count_; ++qi) {
        const auto& q=qpoint_fractions_[qi];
        const auto k=mode_wavevector({qi,0});
        const double k2=k[0]*k[0]+k[1]*k[1]+k[2]*k[2];
        int lo[3],hi[3]; double candidates=1;
        for (int a=0; a<3; ++a) {
            const double extent=std::sqrt(k2*(inv[a][0]*inv[a][0]+inv[a][1]*inv[a][1]+inv[a][2]*inv[a][2]));
            if (!std::isfinite(extent) || extent>100 || std::abs(q[a])>100)
                throw std::invalid_argument("Callaway first-BZ validation exceeds supported lattice conditioning.");
            lo[a]=static_cast<int>(std::ceil(q[a]-extent-1e-10));
            hi[a]=static_cast<int>(std::floor(q[a]+extent+1e-10));
            candidates*=hi[a]-lo[a]+1;
        }
        if (candidates>100000) throw std::invalid_argument("Callaway first-BZ validation work limit exceeded.");
        for (int i=lo[0]; i<=hi[0]; ++i) for (int j=lo[1]; j<=hi[1]; ++j) for (int l=lo[2]; l<=hi[2]; ++l) {
            const double shift[3]={q[0]-i,q[1]-j,q[2]-l};
            double length2=0;
            for (int a=0; a<3; ++a) {
                const double value=b[a][0]*shift[0]+b[a][1]*shift[1]+b[a][2]*shift[2];
                length2+=value*value;
            }
            if (length2<k2-1e-5*std::max(k2,1e-20))
                throw std::invalid_argument("Callaway qpoint is outside the first Brillouin zone; supply first-BZ representatives.");
        }
    }
    size_t hi=static_cast<size_t>(std::lower_bound(temperature_samples_.begin(),temperature_samples_.end(),temperature)-temperature_samples_.begin());
    const size_t lo=hi>0 ? hi-1 : hi;
    const double fraction=hi==lo ? 0 : (temperature-temperature_samples_[lo])/(temperature_samples_[hi]-temperature_samples_[lo]);
    const size_t nm=static_cast<size_t>(qpoint_count_)*branch_count_;
    phonomc::CallawayData data;
    size_t positive=0;
    for (double omega : mode_angular_frequency_data_) if (omega>0) ++positive;
    if (positive!=active_mode_list_.size())
        throw std::invalid_argument("Callaway requires all positive-frequency modes, including stationary modes.");
    for (const auto& mode : active_mode_list_) {
        const size_t flat=static_cast<size_t>(flatten_mode_index(mode));
        data.omega.push_back(mode_angular_frequency(mode));
        data.capacity.push_back(mode_heat_capacity(temperature,mode));
        data.wavevector.push_back(mode_wavevector(mode));
        const auto rate=[&](const std::vector<double>& table) {
            return 4*std::acos(-1.0)*((1-fraction)*table[lo*nm+flat]+fraction*table[hi*nm+flat]);
        };
        data.normal_rate.push_back(rate(gamma_normal_));
        data.resistive_rate.push_back(rate(gamma_umklapp_) + rate(gamma_extra_));
    }
    return data;
}

std::array<double,2> PhononMaterial::mode_callaway_rates(double temperature,const Mode& mode) const {
    if(gamma_normal_.empty() || gamma_umklapp_.empty() || !std::isfinite(temperature) || temperature<=0 ||
       temperature<temperature_samples_.front() || temperature>temperature_samples_.back())
        throw std::invalid_argument("Nonlinear Callaway temperature " + std::to_string(temperature) + " K is outside the supplied N/U table");
    if(mode[0]<0 || mode[0]>=qpoint_count_ || mode[1]<0 || mode[1]>=branch_count_)
        throw std::out_of_range("Invalid Callaway mode");
    size_t hi=std::lower_bound(temperature_samples_.begin(),temperature_samples_.end(),temperature)-temperature_samples_.begin();
    size_t lo=hi?hi-1:hi;double f=lo==hi?0:(temperature-temperature_samples_[lo])/(temperature_samples_[hi]-temperature_samples_[lo]);
    const size_t m=flatten_mode_index(mode),nm=static_cast<size_t>(qpoint_count_)*branch_count_;
    auto rate=[&](const std::vector<double>& a){return 4*std::acos(-1.)*((1-f)*a[lo*nm+m]+f*a[hi*nm+m]);};
    return {rate(gamma_normal_),rate(gamma_umklapp_)+rate(gamma_extra_)};
}
