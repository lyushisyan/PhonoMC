#pragma once

#include "SimulationConfig.h"
#include "SurfaceMesh.h"

#include <array>
#include <string>
#include <vector>

class SimulationDomain {
public:
    explicit SimulationDomain(const SimulationConfig& args);

    const SurfaceMesh& mesh() const { return mesh_; }
    int subvolume_count() const { return subvolume_count_; }
    double volume() const { return volume_; }
    const std::array<double, 3>& bounds_min() const { return bounds_min_; }
    const std::array<double, 3>& bounds_max() const { return bounds_max_; }
    const std::vector<std::array<double, 3>>& subvolume_centers() const { return subvolume_centers_; }
    const std::vector<double>& subvolume_volumes() const { return subvolume_volumes_; }
    const std::vector<std::array<int, 2>>& subvolume_connections() const { return subvolume_connections_; }
    const std::vector<char>& boundary_conditions() const { return facet_boundary_conditions_; }
    const std::vector<int>& reservoir_facets() const { return reservoir_facets_; }
    const std::vector<double>& reservoir_values() const { return reservoir_values_; }
    char facet_boundary_condition(int facet) const;
    double reservoir_value_for_facet(int facet, double fallback = 300.0) const;
    bool has_periodic_pair(int facet) const;
    int periodic_pair_facet(int facet) const;
    std::array<double, 3> periodic_shift_for_facet(int facet) const;
    bool is_rough_facet(int facet) const;
    double roughness_for_facet(int facet, double fallback = 0.0) const;

private:
    void build_surface_mesh(const SimulationConfig& args);
    void sync_surface_mesh_properties();
    void assign_boundary_conditions(const SimulationConfig& args);
    void build_periodic_connections(const SimulationConfig& args);
    void initialize_subvolumes(const SimulationConfig& args);
    void initialize_slice_subvolumes(const SimulationConfig& args);
    void initialize_grid_subvolumes(const SimulationConfig& args);
    void build_subvolume_connections();
    void write_domain_summary(const SimulationConfig& args) const;

    SurfaceMesh mesh_;
    std::array<double, 3> bounds_min_ {0.0, 0.0, 0.0};
    std::array<double, 3> bounds_max_ {1.0, 1.0, 1.0};

    int facet_count_ = 0;
    int subvolume_count_ = 1;
    double volume_ = 1.0;

    std::vector<char> facet_boundary_conditions_;
    std::vector<int> boundary_facets_;
    std::vector<int> reservoir_facets_;
    std::vector<double> reservoir_values_;
    std::vector<int> rough_facets_;
    std::vector<double> roughness_values_;
    std::vector<std::array<int, 2>> connected_facets_;
    std::vector<std::array<double, 3>> subvolume_centers_;
    std::vector<double> subvolume_volumes_;
    std::vector<std::array<int, 2>> subvolume_connections_;
    std::vector<int> periodic_pair_;
    std::vector<std::array<double, 3>> periodic_shift_;
};
