#include "solver/ChannelCell.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace phonomc::channel {
namespace {
Vec sub(Vec a,Vec b){for(int j=0;j<3;++j)a[j]-=b[j];return a;}
double dot(Vec a,Vec b){double v=0;for(int j=0;j<3;++j)v+=a[j]*b[j];return v;}
Vec cross(Vec a,Vec b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
Vec mean(const Face& f){Vec c{};for(auto p:f)for(int j=0;j<3;++j)c[j]+=p[j]/f.size();return c;}
bool close(Vec a,Vec b){return dot(sub(a,b),sub(a,b))<1e-24;}
void append(Face& f,Vec p){if(f.empty()||!close(f.back(),p))f.push_back(p);}
Polyhedron clip(Polyhedron poly,Vec normal,Vec origin,double sign){
    Polyhedron result;Face cap;
    const auto distance=[&](Vec p){return sign*dot(normal,sub(p,origin))-.5;};
    for(const auto& face:poly){
        Face out;
        for(size_t i=0;i<face.size();++i){
            Vec a=face[i],b=face[(i+1)%face.size()];double da=distance(a),db=distance(b);
            const bool ia=da<=1e-12,ib=db<=1e-12;
            if(ia)append(out,a);
            if(ia!=ib){
                double t=da/(da-db);Vec x{};for(int j=0;j<3;++j)x[j]=a[j]+t*(b[j]-a[j]);
                append(out,x);if(std::none_of(cap.begin(),cap.end(),[&](Vec c){return close(c,x);}))cap.push_back(x);
            }
        }
        if(out.size()>1&&close(out.front(),out.back()))out.pop_back();
        if(out.size()>=3)result.push_back(std::move(out));
    }
    if(cap.size()>=3){
        Vec c=mean(cap),n=normal;double norm=std::sqrt(dot(n,n));for(double& v:n)v/=norm;
        Vec u=cross(n,std::abs(n[0])<.8?Vec{1,0,0}:Vec{0,1,0});norm=std::sqrt(dot(u,u));for(double& v:u)v/=norm;
        Vec v=cross(n,u);
        std::sort(cap.begin(),cap.end(),[&](Vec a,Vec b){a=sub(a,c);b=sub(b,c);return std::atan2(dot(a,v),dot(a,u))<std::atan2(dot(b,v),dot(b,u));});
        result.push_back(std::move(cap));
    }
    return result;
}
}
Cell::Cell(Vec origin,Mat b):center(origin),basis(b){
    Vec c0{b[0][0],b[1][0],b[2][0]},c1{b[0][1],b[1][1],b[2][1]},c2{b[0][2],b[1][2],b[2][2]};
    double det=dot(c0,cross(c1,c2));volume=std::abs(det);
    if(!std::isfinite(volume)||volume<=0)throw std::invalid_argument("Degenerate phonon channel cell");
    inverse={cross(c1,c2),cross(c2,c0),cross(c0,c1)};for(auto& row:inverse)for(double& x:row)x/=det;
    for(int i=0;i<3;++i){double radius=0;for(int j=0;j<3;++j)radius+=.5*std::abs(b[i][j]);lower[i]=center[i]-radius;upper[i]=center[i]+radius;}
}
Vec Cell::coordinates(Vec point)const{Vec y{};for(int i=0;i<3;++i)y[i]=dot(inverse[i],sub(point,center));return y;}
Polyhedron Cell::polyhedron()const{
    std::array<Vec,8> vertices;for(int i=0;i<8;++i){vertices[i]=center;for(int d=0;d<3;++d)for(int j=0;j<3;++j)vertices[i][d]+=basis[d][j]*((i&(1<<j))?.5:-.5);}
    constexpr int faces[6][4]={{0,2,6,4},{1,3,7,5},{0,1,5,4},{2,3,7,6},{0,1,3,2},{4,5,7,6}};
    Polyhedron poly;for(const auto& f:faces){Face face;for(int i:f)face.push_back(vertices[i]);poly.push_back(std::move(face));}return poly;
}
Polyhedron intersect(const Cell& a,const Cell& b){
    for(int i=0;i<3;++i)if(a.upper[i]<=b.lower[i]||b.upper[i]<=a.lower[i])return {};
    auto poly=a.polyhedron();for(int i=0;i<3&&!poly.empty();++i)for(double s:{-1.,1.})poly=clip(std::move(poly),b.inverse[i],b.center,s);
    return poly;
}
std::vector<QuadraturePoint> quadrature(const Polyhedron& poly){
    Face all;for(const auto& f:poly)for(auto x:f)all.push_back(x);if(all.empty())return {};
    const Vec anchor=mean(all);std::vector<QuadraturePoint> points;
    constexpr double large=.5854101966249685,small=.1381966011250105;
    for(const auto& face:poly)for(size_t i=1;i+1<face.size();++i){
        std::array<Vec,4> p{anchor,face[0],face[i],face[i+1]};
        double vol=std::abs(dot(sub(p[1],p[0]),cross(sub(p[2],p[0]),sub(p[3],p[0]))))/6;
        if(vol<=0)continue;
        for(int a=0;a<4;++a){Vec x{};for(int b=0;b<4;++b)for(int j=0;j<3;++j)x[j]+=(a==b?large:small)*p[b][j];points.push_back({x,vol/4});}
    }
    return points;
}
}
