#include "solver/AcousticInterfaceSampler.h"
#include "PhononMaterial.h"
#include "solver/ChannelCell.h"
#include <memory>
#include <limits>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <string>

namespace phonomc {
double AcousticInterfaceSampler::transmission(double z1,double z2,double c1,double c2) {
    if (!std::isfinite(z1)||!std::isfinite(z2)||!std::isfinite(c1)||!std::isfinite(c2)||
        z1<=0||z2<=0||c1<0||c2<0||c1>1||c2>1)
        throw std::invalid_argument("Invalid scalar acoustic impedance/angle.");
    const double a=z1*c1,b=z2*c2;
    return a+b>0 ? std::clamp(4*a*b/((a+b)*(a+b)),0.0,1.0) : 0;
}

void AcousticInterfaceSampler::configure(const std::vector<const PhononMaterial*>& mats,int axis,
    const std::vector<double>& densities,double df,double dk) {
    if (axis<0||axis>2||mats.size()<2||densities.size()!=mats.size()||
        !std::isfinite(df)||df<=0||!std::isfinite(dk)||dk<=0)
        throw std::invalid_argument("AMM needs layered materials, mass densities and positive frequency/k-parallel bins.");
    const double pi=std::acos(-1.0);
    const int a=(axis+1)%3,b=(axis+2)%3;
    struct Mode { double omega,flux; PhononMaterial::Vec3 k,v; int branch; };
    std::vector<std::vector<Mode>> modes(mats.size());
    std::vector<std::vector<int>> mirror(mats.size());
    for (size_t m=0;m<mats.size();++m) {
        if (!mats[m]||!std::isfinite(densities[m])||densities[m]<=0)
            throw std::invalid_argument("AMM needs finite positive material densities.");
        const auto& mat=*mats[m];const auto& mesh=mat.q_mesh();
        if(mesh[0]<=0||mesh[1]<=0||mesh[2]<=0)
            throw std::runtime_error("AMM needs a complete uniform q mesh.");
        const channel::Cell reciprocal({0,0,0},mat.reciprocal_lattice());
        using GridKey=std::array<long long,3>;
        const auto grid_key=[&](PhononMaterial::Vec3 k){
            auto q=reciprocal.coordinates(k);GridKey key{};
            for(int j=0;j<3;++j){
                const double x=q[j]*2*mesh[j];const long long index=std::llround(x),period=2LL*mesh[j];
                if(std::abs(x-index)>1e-5)throw std::runtime_error("Interface reflection does not preserve the material q grid.");
                key[j]=(index%period+period)%period;
            }
            return key;
        };
        const int n=mat.active_mode_count();mirror[m].assign(n,-1);
        std::map<GridKey,int> qpoints;
        for(int i=0;i<n;++i){
            const auto id=mat.active_mode_at(i);const auto v=mat.mode_group_velocity(id);
            const double w=mat.mode_angular_frequency(id);const auto k=mat.mode_wavevector(id);
            modes[m].push_back({w,std::abs(v[axis])/mat.energy_density_normalization(),k,v,id[1]});
            qpoints[grid_key(k)]=id[0];
        }
        std::vector<int> optical_negative;
        for(int j=0;j<n;++j)if(modes[m][j].branch>=3 && modes[m][j].v[axis]<-1e-7)optical_negative.push_back(j);
        std::sort(optical_negative.begin(),optical_negative.end(),[&](int i,int j){return std::tie(modes[m][i].omega,i)<std::tie(modes[m][j].omega,j);});
        // Locate the reflected reciprocal-grid address, then resolve degenerate
        // branches by frequency and velocity. Rounded frequency/velocity buckets
        // are not used: a bucket boundary must not break a valid symmetry pair.
        for(int i=0;i<n;++i){
            const auto& x=modes[m][i];if(x.v[axis]<=1e-7)continue;
            auto reflected_k=x.k;reflected_k[axis]*=-1;
            auto it=qpoints.find(grid_key(reflected_k));
            if(it==qpoints.end())throw std::runtime_error("AMM reflection has no reciprocal-grid partner.");
            int best=-1;double error=std::numeric_limits<double>::infinity();
            for(int branch=0;branch<mat.branch_count();++branch){
                if((x.branch<3)!=(branch<3))continue;
                const int j=mat.active_index_for_mode({it->second,branch});if(j<0||mirror[m][j]>=0)continue;
                const auto& y=modes[m][j];if(y.v[axis]>=-1e-7)continue;
                if(std::abs(x.omega-y.omega)>1e-6*std::max(x.omega,y.omega))continue;
                double dv=0,scale=0;for(int d=0;d<3;++d){dv+=std::pow(x.v[d]-(d==axis?-y.v[d]:y.v[d]),2);scale+=x.v[d]*x.v[d];}
                if(dv>1e-12*std::max(1.,scale))continue;
                // Prefer the original branch when degenerate velocities tie.
                dv+=(branch==x.branch?0:1e-24);
                if(dv<error){error=dv;best=j;}
            }
            // Optical states never transmit in this scalar acoustic model.
            // Some nearly degenerate optical eigenvectors are not symmetry
            // covariant at one q point; find an actual frequency/velocity mirror
            // elsewhere in that optical branch, with the same strict tolerance.
            if(best<0 && x.branch>=3){
                const double tol=1e-6*x.omega;
                auto first=std::lower_bound(optical_negative.begin(),optical_negative.end(),x.omega-tol,
                    [&](int j,double w){return modes[m][j].omega<w;});
                for(auto jt=first;jt!=optical_negative.end() && modes[m][*jt].omega<=x.omega+tol;++jt){
                    const int j=*jt;const auto& y=modes[m][j];if(mirror[m][j]>=0||y.branch!=x.branch)continue;
                    double dv=0,scale=0;for(int d=0;d<3;++d){dv+=std::pow(x.v[d]-(d==axis?-y.v[d]:y.v[d]),2);scale+=x.v[d]*x.v[d];}
                    if(dv<=1e-12*std::max(1.,scale)&&dv<error){error=dv;best=j;}
                }
            }
            if(best<0)throw std::runtime_error("AMM needs mirror-symmetric frequency/velocity data for the interface axis (material "+std::to_string(m)+", active "+std::to_string(i)+", v "+std::to_string(x.v[0])+","+std::to_string(x.v[1])+","+std::to_string(x.v[2])+").");
            mirror[m][i]=best;mirror[m][best]=i;
        }
        for(int i=0;i<n;++i)if(std::abs(modes[m][i].v[axis])>=1e-7&&mirror[m][i]<0)
            throw std::runtime_error("AMM reflection pairing is incomplete.");
    }
    // A q-grid point represents a finite reciprocal-space cell. Linearize
    // omega(k) in that cell and map it to (frequency, k_parallel). Geometric
    // overlaps, not coincident point bins, define available elastic channels.
    using channel::Vec;using channel::Mat;using channel::Cell;
    std::vector<Mat> qbasis(mats.size());
    std::vector<std::vector<std::unique_ptr<Cell>>> cells(mats.size());
    for(size_t m=0;m<mats.size();++m) {
        const auto& mesh=mats[m]->q_mesh();const auto& rec=mats[m]->reciprocal_lattice();
        if(mesh[0]<=0||mesh[1]<=0||mesh[2]<=0||
           static_cast<long long>(mesh[0])*mesh[1]*mesh[2]!=mats[m]->qpoint_count())
            throw std::runtime_error("Cell-overlap AMM requires a complete uniform q mesh.");
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)qbasis[m][i][j]=rec[i][j]/mesh[j];
        cells[m].resize(modes[m].size());
        for(size_t i=0;i<modes[m].size();++i) {
            const auto& x=modes[m][i];if(x.branch>=3||std::abs(x.v[axis])<1e-7)continue;
            Mat projected{};
            for(int j=0;j<3;++j){
                for(int d=0;d<3;++d)projected[0][j]+=x.v[d]*qbasis[m][d][j]/(2*pi);
                projected[1][j]=qbasis[m][a][j];projected[2][j]=qbasis[m][b][j];
            }
            cells[m][i]=std::make_unique<Cell>(Vec{x.omega/(2*pi),x.k[a],x.k[b]},projected);
        }
    }
    struct Key {int branch,f,x,y;bool operator==(const Key& k)const{return std::tie(branch,f,x,y)==std::tie(k.branch,k.f,k.x,k.y);}};
    struct Hash {size_t operator()(const Key& k)const{size_t h=0;for(int v:{k.branch,k.f,k.x,k.y})h^=std::hash<int>{}(v)+0x9e3779b9+(h<<6)+(h>>2);return h;}};
    const Vec bin{df,dk,dk};
    const auto bounds=[&](const Cell& c){
        std::array<std::array<int,2>,3> bs{};long double count=1;
        for(int j=0;j<3;++j){
            const double lo=std::floor(c.lower[j]/bin[j]),hi=std::floor(c.upper[j]/bin[j]);
            if(!std::isfinite(lo)||!std::isfinite(hi)||lo<std::numeric_limits<int>::min()||hi>=std::numeric_limits<int>::max())
                throw std::runtime_error("AMM search grid is too fine; enlarge search bins (this does not change channel integrals).");
            bs[j]={static_cast<int>(lo),static_cast<int>(hi)};count*=hi-lo+1;
        }
        if(count>1000000)throw std::runtime_error("AMM search bins create too many candidates per cell; enlarge search bins.");
        return bs;
    };
    std::vector<Interface> built(mats.size()-1);
    for(size_t face=0;face+1<mats.size();++face) {
        auto& dst=built[face];std::array<std::vector<int>,2> incoming;
        std::array<std::vector<double>,2> load;
        for(int side=0;side<2;++side){
            const size_t m=face+side;dst.rows[side].resize(modes[m].size());load[side].resize(modes[m].size());
            for(size_t i=0;i<modes[m].size();++i){const auto& x=modes[m][i];
                if(side==0?x.v[axis]<=1e-7:x.v[axis]>=-1e-7)continue;
                auto& r=dst.rows[side][i];r.reflected=mirror[m][i];r.flux=x.flux;
                if(r.reflected<0)throw std::runtime_error("Missing AMM reflection partner.");
                if(cells[m][i]){incoming[side].push_back(i);dst.audit.incident_acoustic_flux[side]+=x.flux;}
            }
        }
        // Identical discrete banks must have no artificial interface in the
        // acoustic subspace, independent of search-bin widths.
        bool identical=modes[face].size()==modes[face+1].size() &&
            mats[face]->energy_density_normalization()==mats[face+1]->energy_density_normalization() &&
            densities[face]==densities[face+1];
        if(identical)for(size_t i=0;i<modes[face].size();++i){const auto& x=modes[face][i];const auto& y=modes[face+1][i];
            if(x.omega!=y.omega||x.k!=y.k||x.v!=y.v){identical=false;break;}}
        if(identical){
            for(int side=0;side<2;++side)for(int i:incoming[side]){
                // Same active state has the outgoing velocity in the other bank.
                dst.rows[side][i].transmitted.push_back({i,1});
                dst.audit.transmitted_flux[side]+=dst.rows[side][i].flux;
            }
            continue;
        }
        const bool direct=incoming[0].size()*incoming[1].size()<=4096;
        std::unordered_map<Key,std::vector<int>,Hash> index;
        if(!direct)for(int j:incoming[1]){const auto& c=*cells[face+1][j];const auto bs=bounds(c);
            for(int f=bs[0][0];f<=bs[0][1];++f)for(int x=bs[1][0];x<=bs[1][1];++x)for(int y=bs[2][0];y<=bs[2][1];++y)
                index[{modes[face+1][j].branch,f,x,y}].push_back(j);
        }
        struct Pair {int left,right;double g;};std::vector<Pair> pairs;
        std::vector<int> seen(modes[face+1].size(),-1);
        const auto impedance=[&](size_t m,int i,Vec point){
            const auto& mode=modes[m][i];const Vec local=cells[m][i]->coordinates(point);Vec k=mode.k;
            for(int d=0;d<3;++d)for(int j=0;j<3;++j)k[d]+=qbasis[m][d][j]*local[j];
            const double norm=std::sqrt(std::inner_product(k.begin(),k.end(),k.begin(),0.0));
            return std::array<double,2>{norm>0?densities[m]*point[0]*2*pi/norm:0,norm>0?std::abs(k[axis])/norm:0};
        };
        for(int i:incoming[0]){
            const auto& c=*cells[face][i];std::vector<int> candidates;
            if(direct){for(int j:incoming[1])if(modes[face][i].branch==modes[face+1][j].branch)candidates.push_back(j);}
            else {
            const auto bs=bounds(c);
            for(int f=bs[0][0];f<=bs[0][1];++f)for(int x=bs[1][0];x<=bs[1][1];++x)for(int y=bs[2][0];y<=bs[2][1];++y){
                auto it=index.find({modes[face][i].branch,f,x,y});if(it==index.end())continue;
                for(int j:it->second)if(seen[j]!=i){seen[j]=i;candidates.push_back(j);}
            }
            }
            // Stable accumulation makes search-grid resolution an implementation
            // choice only: it cannot open/close a physical transmission channel.
            std::sort(candidates.begin(),candidates.end());
            for(int j:candidates){
                const auto& other=*cells[face+1][j];const auto points=channel::quadrature(channel::intersect(c,other));
                if(points.empty())continue;
                double overlap=0,g=0;
                const double density=std::min(modes[face][i].flux/c.volume,modes[face+1][j].flux/other.volume);
                for(const auto& q:points){
                    if(q.point[0]<=0)continue;
                    auto z1=impedance(face,i,q.point),z2=impedance(face+1,j,q.point);
                    overlap+=q.weight*density;
                    if(z1[0]>0&&z2[0]>0)g+=q.weight*density*transmission(z1[0],z2[0],std::min(1.,z1[1]),std::min(1.,z2[1]));
                }
                dst.audit.geometric_overlap_flux+=overlap;
                if(g>0){pairs.push_back({i,j,g});load[0][i]+=overlap;load[1][j]+=overlap;}
            }
        }
        // Linearized neighboring cells can overlap near curved/turning
        // dispersion surfaces. Bound BOTH marginals with the same edge factor,
        // preserving reciprocity. Record the correction rather than hiding it.
        for(const auto& pair:pairs){
            const double wl=modes[face][pair.left].flux,wr=modes[face+1][pair.right].flux;
            const double factor=std::max({1.,load[0][pair.left]/wl,load[1][pair.right]/wr});
            dst.audit.maximum_capacity_ratio=std::max(dst.audit.maximum_capacity_ratio,factor);
            const double g=pair.g/factor;dst.audit.capacity_removed_flux+=pair.g-g;
            dst.rows[0][pair.left].transmitted.push_back({mirror[face+1][pair.right],g/wl});
            dst.rows[1][pair.right].transmitted.push_back({mirror[face][pair.left],g/wr});
            dst.audit.transmitted_flux[0]+=g;dst.audit.transmitted_flux[1]+=g;
        }
        for(int side=0;side<2;++side)for(int i:incoming[side]){
            const auto& r=dst.rows[side][i];double sum=0;for(const auto& e:r.transmitted)sum+=e.probability;
            if(sum>1+1e-10)throw std::runtime_error("AMM row exceeds unit probability.");
            if(r.transmitted.empty())dst.audit.no_overlap_flux[side]+=r.flux;
        }
    }
    interfaces_=std::move(built);
}
const AcousticInterfaceSampler::Row& AcousticInterfaceSampler::row(int source,int destination,int active) const {
    if(source<0||destination<0||std::abs(source-destination)!=1)
        throw std::invalid_argument("AMM requires adjacent layers.");
    return interfaces_.at(std::min(source,destination)).rows[source<destination?0:1].at(active);
}
AcousticInterfaceSampler::Outcome AcousticInterfaceSampler::sample(int source,int destination,int active,std::mt19937_64& rng) const {
    const auto& r=row(source,destination,active);
    if(r.reflected<0)throw std::runtime_error("AMM incoming state is tangent or unsupported.");
    double x=std::uniform_real_distribution<double>(0,1)(rng);
    for(const auto& e:r.transmitted) {x-=e.probability;if(x<0)return {destination,e.active};}
    return {source,r.reflected};
}
}
