#include "solver/TemperatureGradientDrive.h"
#include "solver/ParticleStorage.h"
#include "SimulationDomain.h"
#include "PhononMaterial.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace phonomc {
void TemperatureGradientDrive::configure(const SimulationConfig& config, const SimulationDomain& geometry,
                                        const PhononMaterial& material, double particle_volume_a3) {
    validate_scattering_config(config);
    if (!config.temperature_gradient || !(particle_volume_a3 > 0) || !std::isfinite(particle_volume_a3))
        throw std::invalid_argument("Invalid temperature-gradient drive configuration.");
    if (!geometry.is_box_geometry() || !geometry.reservoir_facets().empty())
        throw std::invalid_argument("Temperature-gradient driving requires a periodic box without reservoirs.");
    int paired = 0;
    for (int f=0; f<geometry.mesh().facet_count(); ++f) {
        if (std::abs(geometry.mesh().facet_normals()[f][config.transport_axis]) > 1e-8) {
            if (geometry.facet_boundary_condition(f) != 'P' || !geometry.has_periodic_pair(f))
                throw std::invalid_argument("Temperature-gradient direction must have paired periodic boundaries.");
            ++paired;
        }
    }
    if (paired < 2) throw std::invalid_argument("Missing periodic boundaries along the imposed gradient.");
    material_ = &material; reference_ = config.background_temperature;
    const auto& modes = material.active_mode_list();
    rate_.clear(); hw_.clear(); background_.clear(); cell_scale_.clear();
    std::vector<double> capacity;
    long double sum=0, absolute=0, total_c=0;
    for (const auto& mode : modes) {
        const double c = material.mode_heat_capacity(reference_,mode);
        // velocity is Angstrom/ps, gradient is K/m.
        const double rate = -c*material.mode_group_velocity(mode)[config.transport_axis]*1e-10*(*config.temperature_gradient);
        if (!(c>0) || !std::isfinite(rate)) throw std::invalid_argument("Invalid gradient source mode.");
        capacity.push_back(c); rate_.push_back(rate);
        hw_.push_back(6.582119569e-4*material.mode_angular_frequency(mode));
        background_.push_back(material.bose_occupation(reference_,mode));
        sum+=rate; absolute+=std::abs(rate); total_c+=c;
    }
    if (modes.empty() || std::abs(sum)>1e-8L*absolute)
        throw std::invalid_argument("Gradient source has nonzero net energy: supply a balanced full-BZ velocity/capacity bank.");
    // Remove only accepted quadrature/roundoff residual in the energy mode.
    // A materially unbalanced bank is rejected above, never silently corrected.
    for (size_t m=0; m<modes.size(); ++m) rate_[m]-=static_cast<double>(sum*capacity[m]/total_c);
    capacity_ = capacity;
    centers_ = geometry.grid_centers();
    const int counts[3]={config.grid.nx,config.grid.ny,config.grid.nz};
    for (int a=0;a<3;++a) half_width_[a]=(geometry.bounds_max()[a]-geometry.bounds_min()[a])/(2*counts[a]);
    for (double volume : geometry.grid_volumes()) cell_scale_.push_back(volume/(modes.size()*particle_volume_a3));
    energy_weight_ = modes.size()*particle_volume_a3/material.energy_density_normalization();
    if (!(energy_weight_>0) || !std::isfinite(energy_weight_)) throw std::invalid_argument("Invalid gradient source normalization.");
}

SourceDeposit TemperatureGradientDrive::warm_start(ParticleStorage& particles, const CellParticleIndex& cells,
                                                  double age, std::mt19937_64& rng) const {
    if (!material_ || !std::isfinite(age) || age<0) throw std::invalid_argument("Invalid warm-start age.");
    if (age==0) return {};
    auto seed=*this;
    const auto data=material_->linearized_rta_data(reference_);
    long double total=0,capacity=0;
    for (size_t m=0;m<rate_.size();++m) {
        const double r=data.resistive_rate[m];
        const double effective_age=r>0?-std::expm1(-r*age)/r:age;
        seed.rate_[m]*=effective_age;
        total+=seed.rate_[m];capacity+=data.capacity[m];
    }
    // This is an initial guess, not elapsed simulation time or steady film data.
    // Remove its uniform temperature component to start at the chosen T0 energy.
    for (size_t m=0;m<rate_.size();++m) seed.rate_[m]-=static_cast<double>(total*data.capacity[m]/capacity);
    return seed.apply(particles,cells,1,rng);
}

SourceDeposit TemperatureGradientDrive::apply(ParticleStorage& particles, const CellParticleIndex& cells,
                                            double dt, std::mt19937_64& rng) const {
    if (!material_ || !std::isfinite(dt) || dt<0 || !particles.aligned() ||
        !cells.matches(particles,static_cast<int>(centers_.size())))
        throw std::invalid_argument("Invalid temperature-gradient drive step/cell index.");
    SourceDeposit result;
    if (dt==0) return result;
    (void)rng; // Driving adds no carriers and consumes no random samples.
    auto next=particles.occupation;
    long double actual_energy=0;
    for (size_t cell=0;cell<centers_.size();++cell) {
        const int begin=cells.offsets()[cell],end=cells.offsets()[cell+1];
        long double rate_sum=0,capacity_sum=0;
        for(int pos=begin;pos<end;++pos) {
            const int i=cells.indices()[pos],m=material_->active_index_for_mode(particles.modes[i]);
            if(m<0 || particles.material_ids[i]!=0 || !particles.alive[i] ||
               particles.grid_ids[i]!=static_cast<int>(cell) || !std::isfinite(particles.occupation[i]))
                throw std::invalid_argument("Invalid sampled gradient carrier.");
            rate_sum+=rate_[m];capacity_sum+=capacity_[m];
        }
        if (begin==end) continue;
        // Finite Monte Carlo quadrature need not contain exactly balanced +/-q.
        // Remove its energy component to retain a zero-net-heat gradient drive.
        // Each carrier keeps its original phase-space weight (no mode-count cap).
        for(int pos=begin;pos<end;++pos) {
            const int i=cells.indices()[pos],m=material_->active_index_for_mode(particles.modes[i]);
            const long double energy=dt*(rate_[m]-rate_sum*capacity_[m]/capacity_sum);
            next[i]+=static_cast<double>(energy/hw_[m]);
            if(!std::isfinite(next[i])) throw std::overflow_error("Gradient source occupation overflow.");
            actual_energy+=hw_[m]*(static_cast<long double>(next[i])-particles.occupation[i]);
            if(energy!=0) ++result.occupation_updates;
        }
    }
    result.energy_ev=static_cast<double>(actual_energy*energy_weight_);
    if(!std::isfinite(result.energy_ev)) throw std::overflow_error("Gradient source energy overflow.");
    std::copy(next.begin(),next.end(),particles.occupation.begin());
    return result;
}
}
