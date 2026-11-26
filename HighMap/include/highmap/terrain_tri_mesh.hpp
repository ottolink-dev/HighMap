/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once
#include <optional>
#include <vector>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>

#include "highmap/array.hpp"
#include "highmap/spatial/uv_tile.hpp"
#include "highmap/boundary.hpp"
#include "highmap/spatial/bbox2d.hpp"

#include <unordered_map>

namespace hmap
{

struct ComputeMode;
struct VirtualArray;

/**
 * @brief Triangle mesh representation of a terrain surface.
 *
 * Stores vertices, triangulation, connectivity, and provides geometry,
 * interpolation, smoothing, and pathfinding utilities.
 */
class TerrainTriMesh
{
public:
  /** @brief Axis-aligned 2D bounding box. */
  using BoundingBox = BBox2D;

  /** @brief Undirected edge between two vertices. */
  struct Edge
  {
    size_t v0, v1;
    Edge(size_t a, size_t b) : v0(std::min(a, b)), v1(std::max(a, b))
    {
    }
    bool operator==(const Edge &other) const
    {
      return v0 == other.v0 && v1 == other.v1;
    }
  };

  /** @brief Hash functor for Edge. */
  struct EdgeHash
  {
    std::size_t operator()(const Edge &e) const
    {
      return std::hash<size_t>()(e.v0) ^ (std::hash<size_t>()(e.v1) << 1);
    }
  };

  /** @brief Triangle defined by three vertex indices. */
  struct Triangle
  {
    size_t a, b, c;
  };

  /** @brief Neighbor vertex information. */
  struct Neighbor
  {
    size_t index;
    float  distance2d;
  };

  /** @brief Vertex adjacency lists. */
  struct NeighborData
  {
    std::vector<std::vector<Neighbor>> adjacency;

    void clear()
    {
      adjacency.clear();
    }

    void resize(size_t n)
    {
      adjacency.resize(n);
    }
  };

  /** @brief Result of a shortest-path computation. */
  struct ShortestPathResult
  {
    std::vector<float>  distance;
    std::vector<size_t> parent;
  };

  /** @brief Configuration for shortest-path edge costs and heuristics. */
  struct MeshPathCostConfig
  {
    /**
     * @brief Cost metric model used for edge traversal.
     */
    enum class Metric
    {
      /**
       * @brief Pure 2D Euclidean distance in the XY plane.
       *
       * Edge cost is simply the 2D segment length:
       * @f$ \text{cost} = d_{2D} @f$
       */
      EUCLIDEAN_2D,

      /**
       * @brief Symmetric elevation change metric.
       *
       * Treats uphill and downhill elevation changes equally:
       * @f$ \text{cost} = |\Delta z| + w_{elev} \cdot z_{curr} @f$
       * where @f$ w_{elev} @f$ is @ref elevation_weight.
       */
      SYMMETRIC_Z,

      /**
       * @brief Asymmetric hydrological flow metric for streams.
       *
       * Strongly penalizes uphill steps while favoring natural
       * downhill flow:
       * @f[
       *   \text{cost} = (1 - r_{elev}) \cdot d_{2D}
       *               + (\Delta z > 0 \;?\; p_{up} \cdot (\Delta
       * z)^{e_{dist}} : |\Delta z|)
       *               + r_{elev} \cdot \max(0, \Delta z)
       * @f] where @f$ r_{elev} @f$ is @ref elevation_ratio, @f$
       * e_{dist} @f$ is
       * @ref distance_exponent, and @f$ p_{up} @f$ is @ref
       * upward_penalization.
       */
      HYDROLOGICAL
    } metric = Metric::EUCLIDEAN_2D;

    /** @brief Weight for absolute vertex elevation in SYMMETRIC_Z mode. */
    float elevation_weight = 1.0f;

    /** @brief Elevation vs distance balance factor in HYDROLOGICAL mode. */
    float elevation_ratio = 0.5f;

    /** @brief Power exponent applied to upward elevation steps in
       HYDROLOGICAL mode. */
    float distance_exponent = 2.0f;

    /** @brief Penalty multiplier for uphill steps in HYDROLOGICAL mode. */
    float upward_penalization = 100.0f;

    /**
     * @brief Additional penalty favoring paths passing through the center
     * of the domain.
     *
     * Adds @f$ w_{center} \cdot (1 - 2 \cdot d_{border}) @f$ to the edge
     * cost.
     */
    float domain_center_weight = 0.0f;

    /**
     * @brief Optional per-vertex cost multipliers applied to incoming edge
     * costs.
     */
    std::vector<float> vertex_cost_multipliers = {};
  };

public:
  /** @brief Construct an empty mesh. */
  TerrainTriMesh() = default;

  /** @brief Construct a mesh from 3D points. */
  TerrainTriMesh(const std::vector<glm::vec3> &ref_points);

  /** @brief Construct a mesh from coordinate arrays. */
  TerrainTriMesh(const std::vector<float> &x,
                 const std::vector<float> &y,
                 const std::vector<float> &z);

  /** @brief Build the Delaunay triangulation. */
  void triangulate_delaunay();

  /** @brief Compute vertex adjacency. */
  void compute_neighbors();

  /** @brief Compute vertex gradients. */
  void compute_gradients();

  /** @brief Relax vertex positions in the XY plane. */
  void relax_xy(float lambda = 0.5f,
                int   iterations = 1,
                bool  preserve_chull = true);

  /** @brief Relax vertex positions in 3D. */
  void relax_xyz(float lambda = 0.5f,
                 int   iterations = 1,
                 bool  preserve_chull = true);

  /** @brief Apply Taubin smoothing. */
  void relax_xyz_taubin(float lambda = 0.5f,
                        float mu = -0.55f,
                        int   iterations = 1,
                        bool  preserve_chull = true);

  /** @brief Remap elevation values. */
  void remap_z(float vmin = 0.f, float vmax = 1.f);

  /** @brief Limit terrain slopes. */
  void slope_limiter(float max_slope, int iterations = 10, float sigma = 0.1f);

  /** @brief Limit terrain slopes with per-vertex thresholds. */
  void slope_limiter(const std::vector<float> &max_slope,
                     int                       iterations = 10,
                     float                     sigma = 0.1f);

  /** @brief Subdivide mesh triangles. */
  void subdivise();

  /** @brief Compute barycentric coordinates. */
  bool barycentric(const glm::vec2 &p,
                   size_t           i0,
                   size_t           i1,
                   size_t           i2,
                   float           &w0,
                   float           &w1,
                   float           &w2) const;

  /** @brief Find the triangle containing a point. */
  int find_triangle(const glm::vec2 &p,
                    int              start_tri = 0,
                    bool             linear_search = false) const;

  /** @brief Return the neighboring triangle across an edge. */
  int neighbor_triangle(int tri_index, int edge_index) const;

  /** @brief Linearly interpolate elevation. */
  float interpolate_z_linear(const glm::vec2 &p,
                             int             &last_tri,
                             float            fill_value = 0.f) const;

  /** @brief Gradient-enhanced linear interpolation. */
  float interpolate_z_linear_gradient(const glm::vec2 &p,
                                      int             &last_tri,
                                      float            fill_value = 0.f,
                                      float gradient_scaling = 1.f) const;

  /** @brief Nearest-neighbor interpolation. */
  float interpolate_z_nearest(const glm::vec2 &p) const;

  /** @brief Approximate nearest-neighbor interpolation. */
  float interpolate_z_nearest_approx(const glm::vec2 &p,
                                     int             &last_tri,
                                     float            fill_value = 0.f) const;

  /** @brief Return the mesh bounding box. */
  BoundingBox get_bbox() const;

  /** @brief Return the elevation range. */
  glm::vec2 get_range_z() const;

  /** @brief Return reference mesh dimensions. */
  glm::vec3 get_reference_lengths() const;

  /** @brief Compute vertex areas. */
  std::vector<float> get_vertex_areas(bool normalized) const;

  /** @brief Return the projected mesh area. */
  float get_reference_area_xy() const;

  /** @brief Return the mesh surface area. */
  float get_reference_area() const;

  /** @brief Return the number of vertices. */
  size_t size() const;

  /** @brief Compute shortest paths to the convex hull. */
  ShortestPathResult compute_shortest_paths_to_hull(
      bool  use_delta_z = false,
      float elevation_weight = 1.f) const;

  /** @brief Compute shortest paths from multiple sources with a cost config. */
  ShortestPathResult compute_shortest_paths(
      const std::vector<std::pair<size_t, float>> &sources,
      const MeshPathCostConfig                    &config) const;

  /** @brief Compute shortest paths from a single source with a cost config. */
  ShortestPathResult compute_shortest_paths(
      size_t                    start,
      const MeshPathCostConfig &config) const;

  /** @brief Compute shortest paths from a source vertex. */
  ShortestPathResult compute_shortest_paths(
      size_t start,
      bool   use_delta_z = false,
      float  elevation_weight = 1.f,
      float  domain_center_weight = 0.f) const;

  /** @brief Compute the shortest path between two vertices with a cost config.
   */
  std::vector<size_t> shortest_path(size_t                    start,
                                    size_t                    end,
                                    const MeshPathCostConfig &config) const;

  /** @brief Compute the shortest path between two vertices. */
  std::vector<size_t> shortest_path(size_t start,
                                    size_t end,
                                    bool   use_delta_z = false,
                                    float  elevation_weight = 1.f,
                                    float  domain_center_weight = 0.f) const;

  /** @brief Reconstruct a path from start to end from a shortest path result.
   */
  std::vector<size_t> reconstruct_path(size_t                    start,
                                       size_t                    end,
                                       const ShortestPathResult &result) const;

  /** @brief Reconstruct a path to the convex hull. */
  std::vector<size_t> path_to_hull(size_t                    start,
                                   const ShortestPathResult &result) const;

  /** @brief Pick a boundary vertex according to domain boundary and weights. */
  size_t pick_mesh_boundary_vertex(DomainBoundary boundary,
                                   std::uint32_t  seed = 0,
                                   bool           favor_boundary_center = true,
                                   bool           favor_lower_elevation = true,
                                   bool           favor_sinks = true) const;

  /** @brief Return mesh vertices. */
  const std::vector<glm::vec3> &get_points() const;

  /** @brief Return mesh vertices. */
  std::vector<glm::vec3> &get_points();

  /** @brief Return mesh triangles. */
  const std::vector<Triangle> &get_triangles() const;

  /** @brief Return convex hull vertex indices. */
  const std::vector<size_t> &get_convex_hull() const;

  /** @brief Return vertex neighbors. */
  const NeighborData &get_neighbors() const;

  /** @brief Export the mesh as an OBJ file. */
  bool export_obj(const std::string &filepath) const;

  /** @brief Return a mesh summary string. */
  std::string info_string() const;

  /** @brief Print mesh information. */
  void print_info() const;

  /** @brief Rasterize the mesh into an array. */
  Array to_array(const glm::ivec2         &shape,
                 const std::vector<float> &values = {},
                 const BBox2D             &bbox = BBox2D::unit()) const;

  /** @brief Export vertices to CSV. */
  void to_csv(const std::string &fname) const;

private:
  std::vector<glm::vec3> points;
  std::vector<Triangle>  triangles;
  std::vector<size_t>    halfedges;
  std::vector<size_t>    convex_hull;
  NeighborData           neighbors;
  std::vector<glm::vec2> gradients;

private:
  glm::vec2 to_xy(const glm::vec3 &p) const;
};

/** @brief Compute a cubic pulse value for each mesh vertex. */
std::vector<float> cubic_pulse(const TerrainTriMesh &mesh);

/** @brief Generate a terrain mesh from a heightmap. */
TerrainTriMesh generate_terrain_tri_mesh_from_heightmap(const Array &z,
                                                        float        max_error,
                                                        int max_triangles = 0,
                                                        int max_points = 0);

/** @brief Generate a random terrain mesh from a heightmap tile. */
TerrainTriMesh generate_terrain_tri_mesh_from_heightmap_random(
    const Array   &z,
    int            control_points_count,
    std::uint32_t  seed,
    const UvTile  &region = UvTile::unit(),
    ICommunicator *p_comm = nullptr);

} // namespace hmap