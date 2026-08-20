/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file drainage_basin.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for DrainageBasin class and hydrology utilities.
 * @copyright Copyright (c) 2026
 */

#pragma once
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/terrain_tri_mesh.hpp"

#include <unordered_map>

namespace hmap
{

/**
 * @class DrainageBasin
 * @brief Represents a drainage basin network constructed on a 3D terrain mesh.
 *
 * This class handles the construction and analysis of hydrological flow
 * networks (receivers, streams, outlets) on a 3D triangular mesh. It allows
 * simulating
 * river network development, calculating Strahler orders, inverting receiver
 * maps, breaching lakes, and updating elevations based on response times.
 */
class DrainageBasin
{
public:
  /**
   * @brief Construct a new Drainage Basin object.
   * @param xyz_ Input 3D points representing the terrain vertices.
   */
  DrainageBasin(std::vector<glm::vec3> xyz_);

  /**
   * @brief Get the 3D coordinates of the terrain vertices.
   * @return A const reference to the vector of 3D coordinates.
   */
  const std::vector<glm::vec3> &get_xyz() const;

  /**
   * @brief Get the number of vertices in the basin.
   * @return The size as a size_t.
   */
  size_t size() const;

  /**
   * @brief Export the drainage basin data to a CSV file.
   * @param filename Path to the output CSV file.
   */
  void to_csv(const std::string &filename) const;

  // --- Geometry / Mesh ---

  /**
   * @brief Get the underlying terrain tri mesh (const).
   * @return A const reference to the TerrainTriMesh.
   */
  const TerrainTriMesh &get_mesh() const;

  /**
   * @brief Get the underlying terrain tri mesh (non-const).
   * @return A reference to the TerrainTriMesh.
   */
  TerrainTriMesh &get_mesh();

  /**
   * @brief Compute the area of each vertex in the mesh.
   * @return A vector containing the computed vertex areas.
   */
  std::vector<float> compute_vertex_areas() const;

  /**
   * @brief Remap the elevation (z) values of the mesh vertices to a target
   * range.
   * @param zmin Minimum target elevation value.
   * @param zmax Maximum target elevation value.
   */
  void remap(float zmin = 0.f, float zmax = 1.f);

  // --- Flow graph construction ---

  /**
   * @brief Compute the flow receiver for each vertex deterministically.
   */
  void compute_receivers();

  /**
   * @brief Compute the flow receiver for each vertex with noise for stochastic
   * variations.
   * @param seed           Seed for random generation.
   * @param noise_strength Strength of the noise added to elevations during
   *                       receiver computation.
   */
  void compute_receivers(unsigned int seed, float noise_strength = 0.25f);

  /**
   * @brief Update the stream tree stochastically.
   * @param seed           Seed for random generation.
   * @param noise_strength Strength of noise.
   */
  void update_stream_tree(unsigned int seed, float noise_strength);

  /**
   * @brief Update the stream tree deterministically.
   */
  void update_stream_tree();

  /**
   * @brief Update cached traversal orders for upstream/downstream computations.
   */
  void update_traversals();

  /**
   * @brief Get the indices of the outlets in the basin.
   * @return A reference to the vector of outlet indices.
   */
  std::vector<size_t> &get_outlets() const;

  /**
   * @brief Set the outlets of the basin.
   * @param outlet_indices Vector containing the new outlet indices.
   */
  void set_outlets(const std::vector<size_t> &outlet_indices);

  /**
   * @brief Get the receiver index of each vertex.
   * @return A const reference to the vector of receiver indices.
   */
  const std::vector<size_t> &get_receivers() const;

  /**
   * @brief Invert the receiver map to build a map of children/downstream
   * receivers.
   */
  void invert_receiver_map();

  // --- Basin topology utilities ---

  /**
   * @brief Compute whether each vertex is a ridge node (no upstream flow).
   * @return A vector of booleans indicating ridge nodes.
   */
  std::vector<bool> compute_is_ridge_node() const;

  /**
   * @brief Compute the Strahler stream order for each vertex.
   * @return A vector containing the Strahler order of each vertex.
   */
  std::vector<size_t> compute_strahler_order() const;

  /**
   * @brief Find subroots of the flow network.
   * @return A pair containing a vector of subroots indices and a boolean
   * status.
   */
  std::pair<std::vector<size_t>, bool> find_subroots();

  /**
   * @brief Get the main channel paths of the flow network.
   * @return A vector of main channels, each represented as a vector of vertex
   *         indices.
   */
  std::vector<std::vector<size_t>> get_main_channels() const;

  /**
   * @brief Remove lakes by draining local depressions.
   * @param subroot Vector of subroot indices to process.
   */
  void remove_lakes(const std::vector<size_t> &subroot);

  // --- Hydrology computations ---

  /**
   * @brief Compute response times of the basin vertices.
   * @param  area_acc    Accumulated area for each vertex.
   * @param  erodibility Erodibility coefficient for each vertex.
   * @param  m_exp       Erodibility exponent.
   * @return             A vector of response times.
   */
  std::vector<float> compute_response_times(
      const std::vector<float> &area_acc,
      const std::vector<float> &erodibility,
      float                     m_exp) const;

  /**
   * @brief Perform flow breaching to resolve depressions.
   */
  void flow_breach();

  /**
   * @brief Compute the breaching paths for depressions.
   * @return A vector of paths, each represented as a vector of 3D points.
   */
  std::vector<std::vector<glm::vec3>> flow_breach_paths();

  /**
   * @brief Update terrain elevations based on response times and uplift.
   * @param  response_times Vector of vertex response times.
   * @param  uplift_rate    Rate of tectonic uplift.
   * @param  max_slope      Maximum allowed slope for each vertex.
   * @return                The maximum change in elevation.
   */
  float update_elevations(const std::vector<float> &response_times,
                          float                     uplift_rate,
                          const std::vector<float> &max_slope);

  /**
   * @brief Accumulate contributing area down the network by outlet.
   * @param area Input area contribution of each vertex.
   * @param acc  Output accumulated area.
   */
  void accumulate_area_by_outlet(const std::vector<float> &area,
                                 std::vector<float>       &acc) const;

  // --- Traversal helpers ---

  /**
   * @brief Get the cached upstream traversal order from a specific outlet.
   * @param  outlet Index of the outlet.
   * @return        Const reference to the vector of vertex indices in upstream
   *                order.
   */
  const std::vector<size_t> &for_each_upstream(size_t outlet) const;

  /**
   * @brief Get the cached downstream traversal order from a specific outlet.
   * @param  outlet Index of the outlet.
   * @return        A pair of reverse iterators for traversing downstream.
   */
  auto for_each_downstream(size_t outlet) const
  {
    const auto &t = traversals.at(outlet);
    return std::make_pair(t.rbegin(), t.rend());
  }

private:
  // --- Geometry ---

  TerrainTriMesh mesh;

  // --- Flow graph ---

  std::vector<size_t>              receivers;
  std::vector<size_t>              roots; // basin ID
  std::vector<std::vector<size_t>> children;
  std::vector<bool>                outlets_mask;
  mutable std::vector<size_t>      cached_outlets;
  mutable bool                     outlets_dirty = true;
  int                              tick = 0;

  // --- Traversal cache ---

  std::unordered_map<size_t, std::vector<size_t>> traversals;

  // --- Constants ---

  const size_t invalid_index = size_t(-1);
};

// --- FUNCTIONS

/**
 * @brief Find local minima along the border of a set of 3D coordinates.
 * @param  xyz Vector of 3D points.
 * @param  eps Tolerance for coordinates comparison.
 * @return     A vector of indices of the border minima.
 */
std::vector<size_t> find_border_minima(const std::vector<glm::vec3> &xyz,
                                       float eps = 1e-6f);

/**
 * @brief Find sinks along the border of a terrain tri mesh.
 * @param  mesh Input terrain triangular mesh.
 * @param  eps  Tolerance for border coordinate comparison.
 * @return      A vector of indices of the border sinks.
 */
std::vector<size_t> find_border_sinks(TerrainTriMesh &mesh, float eps = 1e-6f);

/**
 * @brief Performs retopology of a heightmap to generate a 3D point cloud.
 * @param  z             Input 2D heightmap array.
 * @param  max_error     Maximum permitted elevation error.
 * @param  max_triangles Maximum number of triangles to generate.
 * @param  max_points    Maximum number of points to generate.
 * @return               A vector of 3D coordinates representing the
 *                       retopologized vertices.
 */
std::vector<glm::vec3> heightmap_retopology(const Array &z,
                                            float        max_error,
                                            int          max_triangles = 0,
                                            int          max_points = 0);

/**
 * @brief Sample a specified number of points along the border.
 * @param  xyz Vector of 3D points.
 * @param  nb  Number of border points to sample.
 * @return     A vector of indices of the sampled border points.
 */
std::vector<size_t> sample_border_points(const std::vector<glm::vec3> &xyz,
                                         size_t                        nb);

} // namespace hmap
