// Geometry/boundary/source preflight only. No particles, collisions or transport.
// Compile with SimulationConfig.cpp, SurfaceMesh.cpp, SimulationDomain.cpp and
// solver/PrescribedHeatSource.cpp; see the model README for the command.
#include "SimulationConfig.h"
#include "SimulationDomain.h"
#include "solver/PrescribedHeatSource.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

void require(bool value, const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
bool near(double a,double b) {return std::abs(a-b)<1e-7*std::max(1.0,std::abs(b));}

int main(int argc,char** argv) {
    try {
        require(argc==6,"Usage: validate_three_layer_fets OUTPUT_JSON FOUR_INPUTS");
        const std::filesystem::path output=argv[1];
        std::filesystem::create_directories(output.parent_path());
        std::ofstream report(output);
        report<<std::setprecision(15)<<"{\"transport_simulated\":false,\"cases\":[\n";
        double reference_power=-1;
        for(int index=2;index<argc;++index) {
            const auto input=std::filesystem::absolute(argv[index]);
            auto config=load_simulation_config(input.string());
            config.output_folder.clear(); // Audit only; do not create simulation-like outputs.
            require(!config.merge_coplanar_facets && config.conservative_boundaries,"Boundary policy differs");
            const SimulationDomain domain(config);
            const bool gaa=input.filename().string().find("gaafet")!=std::string::npos;
            require(domain.grid_count()==86400,"Unexpected occupied-cell count");
            require(near(domain.volume()/1000,43200),"Unexpected physical volume");
            double bottom=0,top=0,sides=0,rough=0;
            int periodic_facets=0;
            const auto& mesh=domain.mesh();
            for(int f=0;f<mesh.facet_count();++f) {
                const auto c=mesh.facet_centroids()[f];
                const auto n=mesh.facet_normals()[f];
                const double area=mesh.facet_areas()[f]/100;
                const char bc=domain.facet_boundary_condition(f);
                const bool is_bottom=near(c[2],0),is_top=near(c[2],500);
                const bool is_side=c[2]<200 &&
                    (near(c[0],0)||near(c[0],520)||near(c[1],0)||near(c[1],240)) && std::abs(n[2])<.1;
                if(is_bottom||is_top) {
                    require(bc=='T'&&near(domain.reservoir_value_for_facet(f),300),"Wrong thermal boundary");
                    if(is_top) {
                        require(c[0]<100 || c[0]>420,"Channel accidentally thermalized");
                        top+=area;
                    } else bottom+=area;
                } else if(is_side) {
                    require(bc=='P'&&domain.has_periodic_pair(f),"Missing substrate periodic pair");
                    const int paired=domain.periodic_pair_facet(f);
                    require(domain.periodic_pair_facet(paired)==f,"Non-reciprocal periodic pair");
                    const auto shift=domain.periodic_shift_for_facet(f);
                    const int axis=std::abs(n[0])>.5?0:1;
                    for(int k=0;k<3;++k)
                        require(near(shift[k],k==axis ? -n[k]*(axis==0?520:240) : 0),"Periodic mapping has spurious tangential shift");
                    // Centroid/edge-spectrum checks alone cannot establish that
                    // an arbitrary surface point maps into the paired triangle.
                    for(int vertex:mesh.faces()[f]) {
                        bool found=false;
                        for(int other:mesh.faces()[paired]) {
                            bool equal=true;
                            for(int k=0;k<3;++k)
                                equal=equal&&near(mesh.vertices()[vertex][k]+shift[k],mesh.vertices()[other][k]);
                            found=found||equal;
                        }
                        require(found,"Periodic triangles do not match under translation");
                    }
                    sides+=area;++periodic_facets;
                } else {
                    require(bc=='R',"Unexpected nonreflecting boundary");
                    require(near(domain.roughness_for_facet(f),10),"Wrong roughness");
                    rough+=area;
                }
            }
            require(near(bottom,1248)&&near(top,480)&&near(sides,3040),"Boundary areas differ");
            phonomc::PrescribedHeatSource source;source.configure(config,domain);
            const double power=source.base_power_w();
            if(reference_power<0)reference_power=power;
            require(std::abs(power-reference_power)<reference_power*1e-11,"Total device powers differ");
            double layer[3]={0,0,0};int heated=0;
            for(int i=0;i<domain.grid_count();++i) {
                require(near(domain.grid_volumes()[i],500),"Cut/incorrect cell volume");
                const double rate=source.cell_rates_evps()[i];
                if(rate==0)continue;
                ++heated;const auto c=domain.grid_centers()[i];
                require(c[0]>360&&c[0]<420,"Heat outside drain-side interval");
                int which=0;
                if(gaa) {
                    require(c[1]>70&&c[1]<170,"Heat outside sheet width");
                    if(c[2]>245&&c[2]<285)which=0;
                    else if(c[2]>330&&c[2]<370)which=1;
                    else if(c[2]>415&&c[2]<455)which=2;
                    else throw std::runtime_error("Heat deposited in gap or outside sheets");
                } else require(c[1]>90&&c[1]<150&&c[2]>200&&c[2]<400,"Heat outside fin");
                layer[which]+=rate/6.241509074e6;
            }
            require(heated==1440,"Unexpected source-cell count");
            for(int j=0;j<(gaa?3:1);++j)
                require(std::abs(layer[j]-power/(gaa?3:1))<power*1e-11,"Unequal layer power");
            if(index>2)report<<",\n";
            report<<"{\"input\":\""<<input.filename().string()<<"\",\"passed\":true,"
                  <<"\"occupied_cells\":"<<domain.grid_count()<<",\"volume_nm3\":"<<domain.volume()/1000
                  <<",\"thermal_bottom_area_nm2\":"<<bottom<<",\"thermal_top_area_nm2\":"<<top
                  <<",\"periodic_substrate_side_area_nm2\":"<<sides
                  <<",\"periodic_facets\":"<<periodic_facets<<",\"period_nm\":[52,24],\"periodic_translation_verified\":true,\"rough_area_nm2\":"<<rough
                  <<",\"source_cells\":"<<heated<<",\"total_power_W\":"<<power<<",\"layer_power_W\":[";
            for(int j=0;j<(gaa?3:1);++j){if(j)report<<",";report<<layer[j];}
            report<<"]}";
        }
        report<<"\n],\"all_passed\":true}\n";
        std::cout<<"Geometry, exact boundary areas, occupied cells and prescribed layer powers passed.\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
