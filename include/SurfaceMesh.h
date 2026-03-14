#pragma once

#include <array>
#include <string>
#include <tuple>
#include <vector>

class SurfaceMesh {
public:
    using Vec3 = std::array<double, 3>;
    using Tri = std::array<int, 3>;

    struct Facet {
        Vec3 normal {};
        Vec3 centroid {};
        double area = 0.0;
        std::vector<int> faces;
    };
    struct FaceQuery {
        int index = -1;
        double distance = 0.0;
        Vec3 closest {};
    };
    struct EdgeQuery {
        int index = -1;
        double distance = 0.0;
        Vec3 closest {};
    };
    struct PointQuery {
        double distance = 0.0;
        Vec3 closest {};
        bool from_edge = false;
    };

    SurfaceMesh() = default;
    SurfaceMesh(std::vector<Vec3> vertices, std::vector<Tri> faces);

    void set_surface_mesh_data(std::vector<Vec3> vertices, std::vector<Tri> faces);
    void shift_to_origin();

    const std::vector<Vec3>& vertices() const { return vertices_; }
    const std::vector<Tri>& faces() const { return faces_; }
    const std::vector<Facet>& facets() const { return facets_; }
    const std::vector<double>& facet_areas() const { return facet_areas_; }
    const std::vector<Vec3>& facet_normals() const { return facet_normals_; }
    const std::vector<Vec3>& facet_centroids() const { return facet_centroids_; }
    const Vec3& bounds_min() const { return bounds_min_; }
    const Vec3& bounds_max() const { return bounds_max_; }
    double volume() const { return volume_; }

    int face_count() const { return static_cast<int>(faces_.size()); }
    int facet_count() const { return static_cast<int>(facets_.size()); }
    int edge_count() const { return static_cast<int>(edges_.size()); }
    int simplex_count() const { return static_cast<int>(simplices_.size()); }
    const std::vector<std::array<int, 2>>& edges() const { return edges_; }
    const std::vector<std::array<int, 2>>& face_adjacency() const { return face_adjacency_; }
    const std::vector<std::array<int, 2>>& facet_adjacency() const { return facets_adjacency_; }
    const std::vector<std::vector<int>>& facet_boundary_edges() const { return facets_boundary_edges_; }
    const std::vector<int>& interfaces() const { return interfaces_; }

    FaceQuery nearest_face(const Vec3& p) const;
    std::vector<FaceQuery> nearest_face(const std::vector<Vec3>& p) const;
    EdgeQuery nearest_edge(const Vec3& p) const;
    std::vector<EdgeQuery> nearest_edge(const std::vector<Vec3>& p) const;
    PointQuery nearest_point(const Vec3& p) const;
    std::vector<PointQuery> nearest_point(const std::vector<Vec3>& p) const;
    int nearest_facet(const Vec3& p) const;
    bool contains_point(const Vec3& p) const;
    std::vector<bool> contains_point(const std::vector<Vec3>& p) const;
    bool contains_point_ray_cast(const Vec3& p) const;
    std::tuple<Vec3, double, int> trace_boundary_intersection(const Vec3& x, const Vec3& v) const;
    std::tuple<std::vector<Vec3>, std::vector<double>, std::vector<int>, std::vector<int>> trace_boundary_intersections(
        const std::vector<Vec3>& x,
        const std::vector<Vec3>& v) const;
    std::vector<Vec3> sample_surface_points(int n, const std::vector<int>& facets = {}) const;
    std::vector<Vec3> sample_volume_points_naive(int n) const;
    std::vector<Vec3> sample_volume_points(int n) const;
    int facet_index_for_face(int face_index) const;
    std::vector<int> facet_indices_for_faces(const std::vector<int>& face_indices) const;
    void build_volume_tetrahedra();
    void erase_vertices(const std::vector<int>& indices);
    void erase_faces(const std::vector<int>& indices);
    void remove_unreferenced_vertices();
    void merge_duplicate_vertices(double tol = 1e-10);
    void write_stl(const std::string& path, const std::string& name) const;

private:
    void rebuild_cached_properties();
    void compute_bounding_box();
    void compute_edge_list();
    void compute_face_metrics();
    void build_face_bvh();
    void compute_face_neighbors();
    void compute_facet_groups();
    void compute_facet_neighbors();
    void compute_interface_faces();
    void compute_enclosed_volume();
    void clear_tetrahedra_cache();
    void ensure_volume_tetrahedra() const;

    struct Simplex {
        std::array<Vec3, 4> v {};
        std::array<Vec3, 2> bounds {};
        std::array<std::array<double, 3>, 3> invA {};
        double volume = 0.0;
    };

    std::vector<Vec3> vertices_;
    std::vector<Tri> faces_;
    std::vector<std::array<int, 2>> edges_;

    std::vector<Vec3> face_normals_;
    std::vector<Vec3> face_centroids_;
    std::vector<double> face_areas_;
    std::vector<double> face_k_;
    std::vector<std::array<Vec3, 2>> face_bounds_;
    struct FaceBvhNode {
        Vec3 bounds_min {};
        Vec3 bounds_max {};
        int left = -1;
        int right = -1;
        int begin = 0;
        int end = 0;
    };
    std::vector<int> face_bvh_indices_;
    std::vector<FaceBvhNode> face_bvh_nodes_;
    int face_bvh_root_ = -1;
    std::vector<int> face_to_facet_;
    std::vector<std::array<int, 3>> face_edges_;
    std::vector<std::vector<int>> edges_faces_;
    std::vector<std::array<int, 2>> face_adjacency_;

    std::vector<Facet> facets_;
    std::vector<double> facet_areas_;
    std::vector<Vec3> facet_normals_;
    std::vector<Vec3> facet_centroids_;
    std::vector<std::vector<int>> facets_edges_;
    std::vector<std::vector<int>> facets_boundary_edges_;
    std::vector<std::array<int, 2>> facets_adjacency_;
    std::vector<int> interfaces_;
    std::vector<int> interfacets_;
    mutable std::vector<Simplex> simplices_;
    mutable bool simplices_ready_ = false;

    Vec3 bounds_min_ {0.0, 0.0, 0.0};
    Vec3 bounds_max_ {1.0, 1.0, 1.0};
    double volume_ = 0.0;
};
