#include "solver/CarrierCollisionOperator.h"
#include "solver/ParticleStorage.h"
#include "PhononMaterial.h"
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>

namespace phonomc {
void CarrierCollisionOperator::configure(const std::vector<const PhononMaterial*>& materials,
                                         double reference, bool rta) {
    if (materials.empty() || !std::isfinite(reference) || reference<=0)
        throw std::invalid_argument("Invalid carrier collision reference/materials.");
    std::vector<CallawayData> tables;
    std::vector<std::vector<double>> background;
    for (const auto* m:materials) {
        if (!m) throw std::invalid_argument("Null carrier collision material.");
        tables.push_back(rta?m->linearized_rta_data(reference):m->callaway_data(reference));
        std::vector<double> bg;
        for (const auto& mode:m->active_mode_list()) bg.push_back(m->bose_occupation(reference,mode));
        background.push_back(std::move(bg));
    }
    materials_=materials;tables_=std::move(tables);background_=std::move(background);
    fallback_cells_.invalidate();
}
double CarrierCollisionOperator::apply(ParticleStorage& p,
    const std::vector<const PhononMaterial*>& materials, const CollisionGridView& grid,
    double dt, double vp) {
    const auto nc=grid.material_ids.size();
    if (materials!=materials_ || tables_.empty() || !std::isfinite(dt) || dt<0 ||
        !std::isfinite(vp) || vp<=0 || !p.aligned() || !nc ||
        grid.particle_counts.size()!=nc || grid.temperatures.size()!=nc ||
        nc>static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Invalid existing-carrier collision input.");
    const CellParticleIndex* cells=grid.cells;
    if (!cells) {fallback_cells_.invalidate();fallback_cells_.ensure(p,static_cast<int>(nc));cells=&fallback_cells_;}
    if (!cells->matches(p,static_cast<int>(nc)) || cells->counts()!=grid.particle_counts)
        throw std::invalid_argument("Carrier collision needs an up-to-date cell index.");
    // Stage all cells. An invalid sample or exponential failure cannot leave a
    // partially updated timestep, including when cells are evaluated in parallel.
    auto next=p.occupation;
    std::vector<long double> residual(nc,0);
    std::vector<std::exception_ptr> errors(nc);
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int cell=0;cell<static_cast<int>(nc);++cell) {
        try {
            const int mid=grid.material_ids[cell];
            if (mid<0 || mid>=static_cast<int>(materials.size()))
                throw std::invalid_argument("Invalid carrier collision material id.");
            const auto& material=*materials[mid];const auto& table=tables_[mid];
            const auto& bg=background_[mid];
            const int begin=cells->offsets()[cell],end=cells->offsets()[cell+1];
            if (begin==end) continue;
            CallawayData sample;std::vector<double> energy;
            const size_t n=end-begin;
            sample.omega.reserve(n);sample.capacity.reserve(n);sample.wavevector.reserve(n);
            sample.normal_rate.reserve(n);sample.resistive_rate.reserve(n);energy.reserve(n);
            std::vector<int> modes;modes.reserve(n);
            for (int pos=begin;pos<end;++pos) {
                const int i=cells->indices()[pos],m=material.active_index_for_mode(p.modes[i]);
                if (m<0 || !p.alive[i] || p.material_ids[i]!=mid || p.grid_ids[i]!=cell ||
                    !std::isfinite(p.occupation[i]))
                    throw std::invalid_argument("Invalid existing collision carrier.");
                modes.push_back(m);sample.omega.push_back(table.omega[m]);
                sample.capacity.push_back(table.capacity[m]);sample.wavevector.push_back(table.wavevector[m]);
                sample.normal_rate.push_back(table.normal_rate[m]);sample.resistive_rate.push_back(table.resistive_rate[m]);
                energy.push_back(6.582119569e-4*table.omega[m]*(p.occupation[i]-bg[m]));
            }
            if (dt==0) continue;
            // Restrict BOTH the collision state and its conservation Gram matrix
            // to this quadrature. Restricting a precomputed full-grid operator
            // would lose energy/momentum through modes absent from the sample.
            const auto evolved=CallawayOperator(sample).evolve(energy,dt);
            for (int pos=begin;pos<end;++pos) {
                const int j=pos-begin,i=cells->indices()[pos],m=modes[j];
                const double hw=6.582119569e-4*table.omega[m];
                next[i]=bg[m]+evolved[j]/hw;
                if (!std::isfinite(next[i])) throw std::overflow_error("Carrier collision occupation overflow.");
                residual[cell]+=hw*(static_cast<long double>(next[i])-p.occupation[i]);
            }
            residual[cell]*=material.active_mode_count()*p.spatial_volume(mid,vp)/material.energy_density_normalization();
        } catch (...) {errors[cell]=std::current_exception();}
    }
    for (const auto& error:errors) if(error) std::rethrow_exception(error);
    long double sum=0;for(auto r:residual)sum+=r;
    if (!std::isfinite(sum)) throw std::overflow_error("Carrier collision energy residual overflow.");
    std::copy(next.begin(),next.end(),p.occupation.begin());
    return static_cast<double>(sum);
}
}
