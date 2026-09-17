#include "solver/NonlinearCallawayCollision.h"
#include "PhononMaterial.h"
#include "solver/ParticleStorage.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
namespace phonomc {
namespace {
constexpr double hbar=6.582119569e-4, kb=8.617333262145e-5;
using Row=std::array<long double,4>;
using Matrix=std::array<Row,4>;
long double bose(long double x) { return x>700 ? 0 : 1/std::expm1(x); }
// Exponential BGK substep: n'=(1-a)n+a*n_target. Match a-weighted
// moments, not rate-weighted moments, to conserve exactly at finite dt.
void channel(const CallawayData& d, std::vector<double>& n, double T,
             const std::vector<double>& rate, double dt, bool normal) {
    const size_t size=n.size();
    std::vector<long double> a(size),metric(size),phi0(size);
    long double activity=0;
    for(size_t i=0;i<size;++i) {
        a[i]=-std::expm1(-rate[i]*dt); activity+=a[i];
        phi0[i]=hbar*d.omega[i]/(kb*T);
        const auto neq=bose(phi0[i]);metric[i]=a[i]*neq*(1+neq);
    }
    if(activity==0) return;
    // Twice weighted-orthogonalize energy and independent momentum constraints.
    // This also handles absent z momentum and sparse/rank-deficient cells.
    std::vector<std::vector<long double>> basis;
    Row lambda{};
    for(int j=0;j<(normal?4:1);++j) {
        std::vector<long double> v(size);
        long double initial=0;
        for(size_t i=0;i<size;++i) {
            v[i]=j==0?phi0[i]:d.wavevector[i][j-1];
            initial+=metric[i]*v[i]*v[i];
        }
        if(initial==0)continue;
        initial=std::sqrt(initial);
        for(auto& x:v)x/=initial;
        for(int pass=0;pass<2;++pass)for(const auto& b:basis) {
            long double dot=0;for(size_t i=0;i<size;++i)dot+=metric[i]*v[i]*b[i];
            for(size_t i=0;i<size;++i)v[i]-=dot*b[i];
        }
        long double length=0;for(size_t i=0;i<size;++i)length+=metric[i]*v[i]*v[i];
        length=std::sqrt(length);if(length<1e-12L)continue;
        for(auto& x:v)x/=length;
        if(j==0)lambda[0]=initial*length;
        basis.push_back(std::move(v));
    }
    const int dim=static_cast<int>(basis.size());
    if(dim==0)throw std::runtime_error("Nonlinear Callaway has no representable energy constraint");
    Row scale{};
    for(int j=0;j<dim;++j)for(size_t i=0;i<size;++i)
        scale[j]+=std::abs(a[i]*basis[j][i])*std::max(static_cast<long double>(n[i]),bose(phi0[i]));
    auto evaluate=[&](const Row& x, Row& g, Matrix* H, std::vector<long double>* target) {
        g={};if(H)*H={};long double objective=0;
        for(size_t i=0;i<size;++i) {
            if(a[i]==0)continue;
            long double phi=0;for(int j=0;j<dim;++j)phi+=x[j]*basis[j][i];
            if(!(phi>0) || !std::isfinite(phi))return std::numeric_limits<long double>::infinity();
            const long double neq=bose(phi);
            if(target)(*target)[i]=neq;
            objective+=a[i]*(-std::log(-std::expm1(-phi))+n[i]*phi);
            for(int j=0;j<dim;++j) {
                g[j]+=a[i]*basis[j][i]*(n[i]-neq);
                if(H)for(int k=0;k<dim;++k)(*H)[j][k]+=a[i]*basis[j][i]*basis[k][i]*neq*(1+neq);
            }
        }
        return objective;
    };
    std::vector<long double> target(size);
    bool converged=false;
    for(int iteration=0;iteration<70;++iteration) {
        Row g{};Matrix H{};const auto objective=evaluate(lambda,g,&H,&target);
        long double error=0;for(int j=0;j<dim;++j)error=std::max(error,std::abs(g[j])/std::max(1e-30L,scale[j]));
        if(error<2e-13L){converged=true;break;}
        Row delta{};for(int j=0;j<dim;++j)delta[j]=-g[j];
        // Small dense Newton solve with partial pivoting.
        for(int j=0;j<dim;++j) {
            int pivot=j;for(int k=j+1;k<dim;++k)if(std::abs(H[k][j])>std::abs(H[pivot][j]))pivot=k;
            if(!(std::abs(H[pivot][j])>1e-30L))throw std::runtime_error("Singular nonlinear Callaway moment Jacobian");
            std::swap(H[j],H[pivot]);std::swap(delta[j],delta[pivot]);
            for(int k=j+1;k<dim;++k){const auto f=H[k][j]/H[j][j];for(int l=j;l<dim;++l)H[k][l]-=f*H[j][l];delta[k]-=f*delta[j];}
        }
        for(int j=dim-1;j>=0;--j){for(int k=j+1;k<dim;++k)delta[j]-=H[j][k]*delta[k];delta[j]/=H[j][j];}
        long double gd=0;for(int j=0;j<dim;++j)gd+=g[j]*delta[j];
        bool accepted=false;
        for(long double step=1;step>1e-14L;step*=.5L) {
            Row trial=lambda,gt{};for(int j=0;j<dim;++j)trial[j]+=step*delta[j];
            const auto value=evaluate(trial,gt,nullptr,nullptr);
            long double next_error=0;for(int j=0;j<dim;++j)next_error=std::max(next_error,std::abs(gt[j])/std::max(1e-30L,scale[j]));
            if(std::isfinite(value) && (value<=objective+1e-4L*step*gd || next_error<error*.9L)) {
                lambda=trial;accepted=true;break;
            }
        }
        if(!accepted)throw std::runtime_error("Nonlinear Callaway moment line search failed");
    }
    if(!converged)throw std::runtime_error("Nonlinear Callaway moment solve did not converge");
    for(size_t i=0;i<size;++i)if(a[i]>0) {
        n[i]=static_cast<double>((1-a[i])*n[i]+a[i]*target[i]);
        if(!std::isfinite(n[i]) || n[i]<0)throw std::runtime_error("Nonlinear Callaway occupation is invalid");
    }
}
}
std::vector<double> nonlinear_callaway_evolve(const CallawayData& d,const std::vector<double>& n,double T,double dt) {
    const size_t size=n.size();
    if(!size || d.omega.size()!=size || d.wavevector.size()!=size || d.normal_rate.size()!=size || d.resistive_rate.size()!=size ||
       !std::isfinite(T)||T<=0||!std::isfinite(dt)||dt<0)throw std::invalid_argument("Invalid nonlinear Callaway state");
    for(size_t i=0;i<size;++i) {
        if(!std::isfinite(n[i])||n[i]<0||!std::isfinite(d.omega[i])||d.omega[i]<=0 ||
           !std::isfinite(d.normal_rate[i])||d.normal_rate[i]<0 || !std::isfinite(d.resistive_rate[i])||d.resistive_rate[i]<0)
            throw std::invalid_argument("Invalid nonlinear Callaway mode");
        for(double q:d.wavevector[i])if(!std::isfinite(q))throw std::invalid_argument("Invalid nonlinear Callaway wavevector");
    }
    auto out=n;if(dt==0)return out;
    channel(d,out,T,d.normal_rate,dt*.5,true);
    channel(d,out,T,d.resistive_rate,dt,false);
    channel(d,out,T,d.normal_rate,dt*.5,true);
    return out;
}
void NonlinearCallawayCollision::configure(const std::vector<const PhononMaterial*>& materials,double reference) {
    if(materials.empty())throw std::invalid_argument("No nonlinear Callaway materials");
    for(const auto* m:materials){if(!m)throw std::invalid_argument("Null material");m->callaway_data(reference);}
    materials_=materials;fallback_cells_.invalidate();
}
double NonlinearCallawayCollision::apply(ParticleStorage& p,const std::vector<const PhononMaterial*>& materials,
    const CollisionGridView& grid,const CollisionStep& step) {
    const size_t nc=grid.material_ids.size();
    if(materials!=materials_ || !p.aligned() || !nc || grid.particle_counts.size()!=nc || grid.temperatures.size()!=nc ||
       !std::isfinite(step.particle_volume_a3)||step.particle_volume_a3<=0)throw std::invalid_argument("Invalid nonlinear collision grid");
    const auto* cells=grid.cells;
    if(!cells){fallback_cells_.invalidate();fallback_cells_.ensure(p,static_cast<int>(nc));cells=&fallback_cells_;}
    if(!cells->matches(p,static_cast<int>(nc)) || cells->counts()!=grid.particle_counts)throw std::invalid_argument("Stale nonlinear collision cell index");
    auto next=p.occupation;std::vector<long double> residual(nc);std::vector<std::exception_ptr> errors(nc);
#ifdef PHONOMC_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for(int cell=0;cell<static_cast<int>(nc);++cell)try {
        const int mid=grid.material_ids[cell];if(mid<0||mid>=static_cast<int>(materials.size()))throw std::invalid_argument("Invalid cell material");
        const auto& m=*materials[mid];const int begin=cells->offsets()[cell],end=cells->offsets()[cell+1];if(begin==end)continue;
        CallawayData d;std::vector<double> n;
        const double T=grid.temperatures[cell],rate_T=step.local_lifetime_temperature?T:step.lifetime_temperature_kelvin;
        for(int pos=begin;pos<end;++pos) {
            int i=cells->indices()[pos];auto mode=p.modes[i];
            if(!p.alive[i]||p.material_ids[i]!=mid||p.grid_ids[i]!=cell)throw std::invalid_argument("Invalid cell carrier");
            auto rates=m.mode_callaway_rates(rate_T,mode);d.omega.push_back(m.mode_angular_frequency(mode));d.wavevector.push_back(m.mode_wavevector(mode));
            d.normal_rate.push_back(rates[0]);d.resistive_rate.push_back(rates[1]);n.push_back(p.occupation[i]);
        }
        auto out=nonlinear_callaway_evolve(d,n,T,step.time_step_ps);
        for(int pos=begin;pos<end;++pos){int j=pos-begin,i=cells->indices()[pos];next[i]=out[j];residual[cell]+=hbar*d.omega[j]*(static_cast<long double>(out[j])-n[j]);}
        residual[cell]*=m.active_mode_count()*p.spatial_volume(mid,step.particle_volume_a3)/m.energy_density_normalization();
    }catch(...){errors[cell]=std::current_exception();}
    for(const auto& e:errors)if(e)std::rethrow_exception(e);
    long double sum=0;for(auto x:residual)sum+=x;
    std::copy(next.begin(),next.end(),p.occupation.begin());return static_cast<double>(sum);
}
}
