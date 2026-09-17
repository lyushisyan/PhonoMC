#include "solver/ModalResampler.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace phonomc {
ResampleResult resample_modal_carriers(ParticleStorage& p, CellParticleIndex& cells,
    int cell_count, const PhononMaterial& material, const BackgroundCache& bg,
    int cap, std::mt19937_64& rng) {
    if (cap < 1 || !p.aligned() || !cells.matches(p,cell_count))
        throw std::invalid_argument("Invalid modal resampling storage/index/cap.");
    const auto& modes=material.active_mode_list();
    const size_t groups=3*modes.size();
    std::vector<int> head(groups,-1),count(groups,0),next(p.size(),-1);
    std::vector<long double> mass(groups,0);
    // Validate before mutation, including energy stored as deviation from fixed T0.
    for (int i=0;i<p.size();++i)
        if (!p.alive[i] || p.material_ids[i]!=0 || !std::isfinite(p.occupation[i]) ||
            material.active_index_for_mode(p.modes[i])<0)
            throw std::invalid_argument("Resampling requires live single-material modal carriers.");
    long double residual=0;
    for (int cell=0;cell<cell_count;++cell) {
        std::fill(head.begin(),head.end(),-1);std::fill(count.begin(),count.end(),0);
        std::fill(mass.begin(),mass.end(),0);
        for (int pos=cells.offsets()[cell];pos<cells.offsets()[cell+1];++pos) {
            int i=cells.indices()[pos],m=material.active_index_for_mode(p.modes[i]);
            const double delta=p.occupation[i]-bg.occupation(p.modes[i]);
            const int sign=delta>0?0:delta<0?1:2;
            const size_t g=3*static_cast<size_t>(m)+sign;
            next[i]=head[g];head[g]=i;++count[g];
            mass[g]+=sign==2?1:std::abs(delta);
        }
        for (size_t g=0;g<groups;++g) {
            if (count[g]<=cap) continue;
            const int sign=static_cast<int>(g%3);
            const auto mode=modes[g/3];const double background=bg.occupation(mode);
            const double hw=6.582119569e-4*material.mode_angular_frequency(mode);
            const long double step=mass[g]/cap;
            const long double offset=std::uniform_real_distribution<double>(0,1)(rng)*step;
            std::vector<std::pair<int,int>> selected;
            long double cumulative=0;int sample=0;
            for (int i=head[g];i>=0;i=next[i]) {
                const double delta=p.occupation[i]-background;
                cumulative+=sign==2?1:std::abs(delta);
                int multiplicity=0;
                while (sample<cap && (offset+sample*step<cumulative || next[i]<0)) {++multiplicity;++sample;}
                if (multiplicity) selected.emplace_back(i,multiplicity);
                p.alive[i]=0;
            }
            long double remaining=sign==2?0:(sign==0?mass[g]:-mass[g]);
            for (size_t j=0;j<selected.size();++j) {
                const auto [i,n]=selected[j];
                const long double target=j+1==selected.size()?remaining:(sign==2?0:(sign==0?1:-1)*step*n);
                p.alive[i]=1;p.occupation[i]=background+static_cast<double>(target);
                const double actual=p.occupation[i]-background;
                p.energies[i]=hw*actual;remaining-=actual;
            }
            // Remaining is old minus new deviation, including final floating-point rounding.
            residual-=hw*remaining*bg.energy_weight();
        }
    }
    ResampleResult result;result.removed=p.compact_alive();
    result.energy_residual_ev=static_cast<double>(residual);
    cells.invalidate();cells.ensure(p,cell_count);
    return result;
}
}
