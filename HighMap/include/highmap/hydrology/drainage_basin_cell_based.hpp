/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file drainage_basin_cell_based.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for DrainageBasinCellBased class.
 * @copyright Copyright (c) 2026
 */

#pragma once
#include <functional>
#include <limits>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @enum FlowDirectionMethod
 * @brief Enumeration of available algorithms for computing flow directions.
 */
enum FlowDirectionMethod : int
{
  FDM_D8,            ///< Standard D8 flow direction algorithm.
  FDM_PRIORITY_FLOOD ///< Priority flood flow routing algorithm.
};

/**
 * @class DrainageBasinCellBased
 * @brief Represents a cell-based hydrology drainage basin network on a 2D
 * heightmap grid.
 *
 * This class implements hydrological flow network computations (receivers,
 * outlets, main channels, upstream traversals) on a regular grid represented by
 * a 2D Array.
 */
class DrainageBasinCellBased
{
public:
  // --- Construction ---

  /**
   * @brief Default constructor.
   */
  DrainageBasinCellBased() = default;

  /**
   * @brief Construct a new cell-based drainage basin using a heightmap.
   * @param z_ The input heightmap array.
   */
  DrainageBasinCellBased(const Array &z_);

  // --- Geometry / Mesh ---

  /**
   * @brief Get the underlying heightmap array.
   * @return A const reference to the heightmap Array.
   */
  const Array &get_z() const;

  // --- Flow graph construction ---

  /**
   * @brief Compute flow receivers using standard D8, with optional noise.
   * @param seed           Seed for random generation.
   * @param noise_strength Strength of elevation noise added for stochastic flow
   *                       paths.
   */
  void compute_receivers(unsigned int seed = 0, float noise_strength = 0.f);

  /**
   * @brief Compute flow receivers using the priority flood routing algorithm.
   */
  void compute_receivers_priority_flood();

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
   * @brief Get the outlets of the basin.
   * @return A vector of 2D grid coordinates of the outlets.
   */
  std::vector<glm::ivec2> get_outlets() const;

  /**
   * @brief Set the outlets of the basin.
   * @param outlet_indices Vector of outlet 2D coordinates.
   */
  void set_outlets(const std::vector<glm::ivec2> &outlet_indices);

  /**
   * @brief Compute upstream traversal orders for the entire grid.
   * @return A vector of paths, each represented as a vector of 2D coordinates.
   */
  std::vector<std::vector<glm::ivec2>> compute_upstream_traversals();

  // --- Basin topology utilities ---

  /**
   * @brief Find subroots of the flow network.
   * @return A pair containing a matrix of subroots and a boolean status.
   */
  std::pair<Mat<glm::ivec2>, bool> find_subroots();

  /**
   * @brief Remove lakes by draining local depressions.
   * @param subroot Matrix representing subroots.
   */
  void remove_lakes(const Mat<glm::ivec2> &subroot);

  /**
   * @brief Get the main channel paths of the flow network.
   * @return A vector of main channels, each represented as a vector of 2D grid
   *         coordinates.
   */
  std::vector<std::vector<glm::ivec2>> get_main_channels() const;

  // --- Hydrology computations ---

  /**
   * @brief Compute response times of the basin cells.
   * @param  area_acc    Accumulated area array.
   * @param  erodibility Erodibility coefficient array.
   * @param  m_exp       Erodibility exponent.
   * @return             An Array of response times.
   */
  Array compute_response_times(const Array &area_acc,
                               const Array &erodibility,
                               float        m_exp) const;

  /**
   * @brief Update elevations based on response times and uplift.
   * @param  response_times Array of cell response times.
   * @param  uplift_rate    Rate of tectonic uplift.
   * @param  max_slope      Maximum allowed slope array.
   * @return                The maximum change in elevation.
   */
  float update_elevations(const Array &response_times,
                          float        uplift_rate,
                          const Array &max_slope);

  /**
   * @brief Accumulate contributing area down the network by outlet.
   * @param acc Output accumulated area array.
   */
  void accumulate_area_by_outlet(Array &acc) const;

  /**
   * @brief Perform flow breaching to resolve depressions.
   */
  void flow_breach();

  // --- Members ---

  Array z; ///< The heightmap array.

  Mat<int>                     outlets_mask; ///< Mask indicating outlet cells.
  Mat<glm::ivec2>              receivers;    ///< Grid of receiver coordinates.
  Mat<std::vector<glm::ivec2>> children;     ///< Grid of children coordinates.
  Mat<glm::ivec2>              roots;        ///< Grid of basin root
                                             // coordinates.

  std::unordered_map<glm::ivec2, std::vector<glm::ivec2>, IVec2Hash>
      traversals; ///<
                  // Cached
                  // traversal
                  // paths.

  const glm::ivec2 null_cell = glm::ivec2(-1, -1); ///< Constant representing an
                                                   // invalid/null cell.

private:
  // constants
  const int   di[8] = {1, 1, 0, -1, -1, -1, 0, 1};
  const int   dj[8] = {0, -1, -1, -1, 0, 1, 1, 1};
  const float cd[8] =
      {1.f, M_SQRT1_2, 1.f, M_SQRT1_2, 1.f, M_SQRT1_2, 1.f, M_SQRT1_2};
};

/**
 * @brief Invert the receiver map to build a map of children/downstream
 * receivers.
 * @param  receivers Input grid of receiver coordinates.
 * @return           A matrix of vectors representing children coordinates for
 *                   each cell.
 */
Mat<std::vector<glm::ivec2>> invert_receiver_map(
    const Mat<glm::ivec2> &receivers);

} // namespace hmap