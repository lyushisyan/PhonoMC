#pragma once
#include <array>
#include <map>
#include <string>
#include <vector>

namespace phonomc {
using MaterialVec3 = std::array<double, 3>;
using MaterialMat3 = std::array<MaterialVec3, 3>;
struct LatticeData {
    double volume_a3 = 0;
    MaterialMat3 reciprocal {};
};
// Owned source data, not a material model. q-major/branch-minor ordering;
// gamma is temperature-major, then q, then branch. No HDF5 handles escape.
struct MaterialData {
    LatticeData lattice;
    int qpoints = 0, branches = 0;
    std::array<int, 3> mesh {};
    std::vector<double> frequency_thz; // raw input, including finite negative values
    std::vector<MaterialVec3> qpoint_fractions, group_velocity; // velocity: Angstrom/ps
    std::vector<double> temperatures, gamma; // K; linewidth in THz (not inverse lifetime)
    // Optional pair, same shape/units as gamma. Required by Callaway.
    std::vector<double> gamma_normal, gamma_umklapp;
    // Independent extra linewidths, never already included in gamma or gamma_U.
    // Each table is (q, branch) or (temperature, q, branch), in THz HWHM.
    std::map<std::string, std::vector<double>> gamma_resistive;
};

// Validate owned raw values at both file and in-memory entry points.
// Throws std::runtime_error; never repairs or changes the input.
void validate_material_data(const MaterialData& data);
} // namespace phonomc
