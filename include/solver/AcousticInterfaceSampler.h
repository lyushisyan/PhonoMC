#pragma once
#include <array>
#include <random>
#include <vector>
class PhononMaterial;
namespace phonomc {
// Dispersive scalar acoustic mismatch approximation. No polarization conversion;
// optical modes reflect. Finite q-cell geometric overlaps and reciprocal flux
// capacities make the discrete kernel stationary. Bin widths only index candidates.
class AcousticInterfaceSampler {
public:
    struct Audit {
        std::array<double,2> incident_acoustic_flux{}, transmitted_flux{}, no_overlap_flux{};
        double geometric_overlap_flux=0, capacity_removed_flux=0, maximum_capacity_ratio=1;
    };
    const Audit& audit(int face) const { return interfaces_.at(face).audit; }
    struct Outcome { int material, active; };
    struct Edge { int active; double probability; };
    struct Row { int reflected=-1; std::vector<Edge> transmitted; double flux=0; };
    void configure(const std::vector<const PhononMaterial*>& materials, int axis,
                   const std::vector<double>& density_kg_m3,
                   double frequency_bin_thz, double parallel_bin_inv_a);
    Outcome sample(int source, int destination, int active, std::mt19937_64& rng) const;
    const Row& row(int source,int destination,int active) const;
    static double transmission(double impedance1,double impedance2,double cosine1,double cosine2);
private:
    struct Interface { std::array<std::vector<Row>,2> rows; Audit audit; };
    std::vector<Interface> interfaces_;
};
}
