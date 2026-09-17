#pragma once
#include <array>
#include <vector>
namespace phonomc::channel {
using Vec = std::array<double,3>;
using Mat = std::array<Vec,3>;
using Face = std::vector<Vec>;
using Polyhedron = std::vector<Face>;
struct Cell {
    Vec center{}, lower{}, upper{};
    Mat basis{}, inverse{};
    double volume=0;
    Cell(Vec origin,Mat full_width_basis);
    Polyhedron polyhedron() const;
    Vec coordinates(Vec point) const;
};
struct QuadraturePoint { Vec point; double weight; };
Polyhedron intersect(const Cell& a,const Cell& b);
// Positive degree-two tetrahedral quadrature; weights are channel-space volumes.
std::vector<QuadraturePoint> quadrature(const Polyhedron& poly);
}
