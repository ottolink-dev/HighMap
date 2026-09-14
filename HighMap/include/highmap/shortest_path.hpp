/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file shortest_path.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Implements shortest path algorithms including Dijkstra's method for 2D
 * array data representation.
 *
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/geometry/path.hpp"

namespace hmap
{

/**
 * @brief Divide the path by adding points based on the lowest elevation
 * difference between each pair of edge endpoints.
 *
 * This method uses the elevation map to subdivide the path by adding
 * intermediate points where the elevation difference between edges is minimal.
 * The `elevation_ratio` parameter balances the influence of absolute elevation
 * versus elevation difference in the cost function used for path finding. The
 * `distance_exponent` affects the weight function used in Dijkstra's algorithm.
 * Areas defined by the `p_mask_nogo` mask are avoided.
 *
 * **Example**
 * @include ex_path_dijkstra.cpp
 *
 * **Result**
 * @image html ex_path_dijkstra.png
 *
 * @param array             Elevation map used to determine elevation
 *                          differences along the path.
 * @param bbox              Bounding box of the domain, defining the area where
 *                          elevation data is valid.
 * @param edge_divisions    Number of subdivisions per edge; set to 0 for
 *                          automatic division based on array shape.
 * @param elevation_ratio   Ratio used to balance absolute elevation and
 *                          elevation difference in the cost function.
 * @param distance_exponent Exponent used in the Dijkstra weight function to
 *                          adjust the influence of distance.
 * @param p_mask_nogo       Optional mask array defining areas to avoid;
 * points in these areas will not be considered.
 *
 * @see                     Array::find_path_dijkstra
 */
Path dijkstra(const Path  &path,
              const Array &array,
              glm::vec4    bbox,
              float        elevation_ratio = 0.f,
              float        distance_exponent = 0.5f,
              float        upward_penalization = 1.f,
              Array       *p_mask_nogo = nullptr);

/**
 * @brief Find a Dijkstra-based cut path between two domain boundaries.
 *
 * Selects the lowest point on the @p start and @p end boundaries of heightmap
 * @p z, then computes a path between them using a weighted Dijkstra search. The
 * path cost combines elevation, distance, and uphill penalization. The result
 * is returned as normalized (x, y) coordinates and elevation samples.
 *
 * @param  z                        Heightmap array.
 * @param  start                    Boundary where the path begins.
 * @param  end                      Boundary where the path ends.
 * @param  dijk_elevation_ratio     Weight of elevation in the cost.
 * @param  dijk_distance_exponent   Exponent applied to distance cost.
 * @param  dijk_upward_penalization Extra penalty for uphill moves.
 *
 * @return                          Path containing normalized (x, y) points and
 *                                  elevations.
 *
 * **Example**
 * @include ex_find_cut_path.cpp
 *
 * **Result**
 * @image html ex_find_cut_path.png
 */
Path find_cut_path_dijkstra(const Array   &z,
                            DomainBoundary start,
                            DomainBoundary end,
                            float          dijk_elevation_ratio = 0.9f,
                            float          dijk_distance_exponent = 2.f,
                            float          dijk_upward_penalization = 100.f,
                            std::uint32_t  seed = 0,
                            bool           favor_boundary_center = true,
                            bool           favor_lower_elevation = true,
                            bool           favor_sinks = true);

/**
 * @brief Generate a stochastic cut path using midpoint displacement.
 *
 * Creates a path between the @p start and @p end boundaries of heightmap @p z
 * by recursively subdividing and perturbing a segment. The result is a natural,
 * irregular path defined by normalized (x, y) coordinates and elevations.
 *
 * @param  z               Heightmap array.
 * @param  start           Start boundary.
 * @param  end             End boundary.
 * @param  seed            Random seed.
 * @param  midp_iterations Number of subdivision steps.
 * @param  midp_sigma      Displacement amplitude.
 *
 * @return                 Path with normalized coordinates and elevations.
 *
 * **Example**
 * @include ex_find_cut_path.cpp
 *
 * **Result**
 * @image html ex_find_cut_path.png
 */
Path find_cut_path_midpoint(const Array   &z,
                            DomainBoundary start,
                            DomainBoundary end,
                            std::uint32_t  seed,
                            float          offset_ratio = 0.2f,
                            int            steps = 16,
                            bool           favor_boundary_center = true,
                            bool           favor_lower_elevation = true,
                            bool           favor_sinks = true);

/**
 * @brief Find a fast coarse-to-fine multiscale cut path between two domain
 * boundaries.
 *
 * Selects candidate boundary cells on the @p start and @p end boundaries, then
 * computes a path using hierarchical multiscale search with corridor
 * restriction.
 *
 * @param  z                        Heightmap array.
 * @param  start                    Boundary where the path begins.
 * @param  end                      Boundary where the path ends.
 * @param  n_levels                 Number of multiscale pyramid levels.
 * @param  corridor_radius          Corridor radius in grid cells.
 * @param  corridor_decay           Corridor radius decay factor per level.
 * @param  dijk_elevation_ratio     Weight of elevation in the cost.
 * @param  dijk_distance_exponent   Exponent applied to distance cost.
 * @param  dijk_upward_penalization Extra penalty for uphill moves.
 * @param  seed                     Seed for boundary cell selection.
 * @param  favor_boundary_center    Bias boundary cell toward center.
 * @param  favor_lower_elevation    Bias boundary cell toward lower elevation.
 * @param  favor_sinks              Bias boundary cell toward local sinks.
 * @param  use_astar                Use A* heuristic.
 * @param  smooth_path              Smooth the resulting path.
 *
 * @return                          Path containing normalized coordinates and
 *                                  sampled elevations.
 */
Path find_cut_path_multiscale(const Array   &z,
                              DomainBoundary start,
                              DomainBoundary end,
                              int            n_levels = 3,
                              int            corridor_radius = 4,
                              float          corridor_decay = 1.f,
                              float          dijk_elevation_ratio = 0.9f,
                              float          dijk_distance_exponent = 2.f,
                              float          dijk_upward_penalization = 100.f,
                              std::uint32_t  seed = 0,
                              bool           favor_boundary_center = true,
                              bool           favor_lower_elevation = true,
                              bool           favor_sinks = true,
                              bool           use_astar = false,
                              bool           smooth_path = false);

/**
 * @brief Finds the path with the lowest elevation and elevation difference
 * between two points in a 2D array using Dijkstra's algorithm.
 *
 * This function calculates the shortest path considering both elevation and
 * elevation differences. It uses a cost function that balances between absolute
 * elevation and elevation change. The path is determined by minimizing the
 * combined cost function.
 *
 * @see                       @cite Dijkstra1971 and
 *                            https://math.stackexchange.com/questions/3088292
 *
 * @warning The `elevation_ratio` parameter must be less than 1 for the
 * algorithm to converge properly.
 *
 * @param ij_start            Starting coordinates (i, j) for the pathfinding.
 * @param ij_end              Ending coordinates (i, j) for the pathfinding.
 * @param i_path[out]         Vector to store the resulting shortest path
 *                            indices in the i direction.
 * @param j_path[out]         Vector to store the resulting shortest path
 *                            indices in the j direction.
 * @param elevation_ratio     Balance factor between absolute elevation and
 *                            elevation difference in the cost function. Should
 *                            be in the range [0, 1[.
 * @param distance_exponent   Exponent used in the distance calculation between
 *                            points. A higher exponent increases the cost of
 *                            elevation gaps, encouraging the path to minimize
 *                            elevation changes and reduce overall cumulative
 *                            elevation gain.
 * @param upward_penalization Penalize upstream slopes.
 * @param p_mask_nogo         Optional pointer to an array mask that defines
 *                            areas to avoid during pathfinding.
 *
 * **Example**
 * @include ex_find_path_dijkstra.cpp
 *
 * **Result**
 * @image html ex_find_path_dijkstra0.png
 * @image html ex_find_path_dijkstra1.png
 * @image html ex_find_path_dijkstra2.png
 */
void find_path_dijkstra(const Array      &z,
                        glm::ivec2        ij_start,
                        glm::ivec2        ij_end,
                        std::vector<int> &i_path,
                        std::vector<int> &j_path,
                        float             elevation_ratio = 0.1f,
                        float             distance_exponent = 2.f,
                        float             upward_penalization = 1.f,
                        const Array      *p_mask_nogo = nullptr);

void find_path_dijkstra(const Array                   &z,
                        glm::ivec2                     ij_start,
                        std::vector<glm::ivec2>        ij_end_list,
                        std::vector<std::vector<int>> &i_path_list,
                        std::vector<std::vector<int>> &j_path_list,
                        float                          elevation_ratio = 0.1f,
                        float                          distance_exponent = 2.f,
                        float        upward_penalization = 1.f,
                        const Array *p_mask_nogo = nullptr);

/**
 * @brief Compute a path between two points using iterative midpoint refinement.
 *
 * The algorithm subdivides segments and shifts midpoints along the
 * perpendicular direction to locally minimize the scalar field `z`. This is a
 * fast heuristic alternative to graph-based methods (not guaranteed optimal).
 *
 * @param  z            2D scalar field (weights).
 * @param  ij_start     Start index.
 * @param  ij_end       End index.
 * @param  offset_ratio Relative transverse displacement (per segment length).
 * @param  max_it       Max iterations (0 = automatic based on distance).
 * @param  steps        Number of samples for transverse search.
 *
 * @return              Vector of grid indices forming the path.
 *
 * **Example**
 * @include ex_find_path_midpoint.cpp
 *
 * **Result**
 * @image html ex_find_path_midpoint.png
 */
std::vector<glm::ivec2> find_path_midpoint(const Array &z,
                                           glm::ivec2   ij_start,
                                           glm::ivec2   ij_end,
                                           float        offset_ratio = 0.2f,
                                           int          max_it = 0,
                                           int          steps = 16);

/**
 * @brief Finds a shortest path on a 2D grid using a fast coarse-to-fine
 * multiscale approach.
 *
 * This function solves pathfinding across multiple grid resolutions to
 * drastically reduce the number of explored cells on large grids while
 * maintaining high path quality. A coarse path is found on downsampled cost
 * maps, upsampled to the next finer resolution, and restricted to a narrow
 * search corridor at each stage.
 *
 * @param z                   Heightmap or cost array.
 * @param ij_start            Starting grid index (i, j).
 * @param ij_end              Target grid index (i, j).
 * @param i_path[out]         Vector storing resulting path i indices.
 * @param j_path[out]         Vector storing resulting path j indices.
 * @param n_levels            Number of multiscale pyramid levels (>= 1).
 * @param corridor_radius     Radius of search corridor in grid cells.
 * @param corridor_decay      Decay factor applied to corridor radius per level.
 * @param elevation_ratio     Weight factor for absolute elevation vs gradient.
 * @param distance_exponent   Exponent applied to elevation differences.
 * @param upward_penalization Penalty factor for moving uphill.
 * @param p_mask_nogo         Optional pointer to mask of impassable regions.
 * @param use_astar           Use A* heuristic instead of Dijkstra.
 * @param smooth_path         Apply smoothing to the resulting path.
 *
 * **Example**
 * @include ex_find_path_multiscale.cpp
 *
 * **Result**
 * @image html ex_find_path_multiscale.png
 */
void find_path_multiscale(const Array      &z,
                          glm::ivec2        ij_start,
                          glm::ivec2        ij_end,
                          std::vector<int> &i_path,
                          std::vector<int> &j_path,
                          int               n_levels = 3,
                          int               corridor_radius = 4,
                          float             corridor_decay = 1.f,
                          float             elevation_ratio = 0.1f,
                          float             distance_exponent = 2.f,
                          float             upward_penalization = 1.f,
                          const Array      *p_mask_nogo = nullptr,
                          bool              use_astar = false,
                          bool              smooth_path = false);

/**
 * @brief Overload of find_path_multiscale returning a vector of grid indices.
 *
 * @param  z                   Heightmap or cost array.
 * @param  ij_start            Starting grid index (i, j).
 * @param  ij_end              Target grid index (i, j).
 * @param  n_levels            Number of multiscale pyramid levels (>= 1).
 * @param  corridor_radius     Radius of search corridor in grid cells.
 * @param  corridor_decay      Decay factor applied to corridor radius per
 *                             level.
 * @param  elevation_ratio     Weight factor for absolute elevation vs gradient.
 * @param  distance_exponent   Exponent applied to elevation differences.
 * @param  upward_penalization Penalty factor for moving uphill.
 * @param  p_mask_nogo         Optional pointer to mask of impassable regions.
 * @param  use_astar           Use A* heuristic instead of Dijkstra.
 * @param  smooth_path         Apply smoothing to the resulting path.
 *
 * @return                     Vector of 2D grid indices forming the path.
 */
std::vector<glm::ivec2> find_path_multiscale(const Array &z,
                                             glm::ivec2   ij_start,
                                             glm::ivec2   ij_end,
                                             int          n_levels = 3,
                                             int          corridor_radius = 4,
                                             float        corridor_decay = 1.f,
                                             float elevation_ratio = 0.1f,
                                             float distance_exponent = 2.f,
                                             float upward_penalization = 1.f,
                                             const Array *p_mask_nogo = nullptr,
                                             bool         use_astar = false,
                                             bool         smooth_path = false);

/**
 * @brief Overload of find_path_multiscale returning a Path object mapped to a
 * bounding box, with point values containing the elevation sampled from the
 * heightmap.
 *
 * @param  z                   Heightmap or cost array.
 * @param  ij_start            Starting grid index (i, j).
 * @param  ij_end              Target grid index (i, j).
 * @param  bbox                Bounding box domain for coordinate remapping.
 * @param  n_levels            Number of multiscale pyramid levels (>= 1).
 * @param  corridor_radius     Radius of search corridor in grid cells.
 * @param  corridor_decay      Decay factor applied to corridor radius per
 *                             level.
 * @param  elevation_ratio     Weight factor for absolute elevation vs gradient.
 * @param  distance_exponent   Exponent applied to elevation differences.
 * @param  upward_penalization Penalty factor for moving uphill.
 * @param  p_mask_nogo         Optional pointer to mask of impassable regions.
 * @param  use_astar           Use A* heuristic instead of Dijkstra.
 * @param  smooth_path         Apply smoothing to the resulting path.
 *
 * @return                     Path containing coordinates and sampled elevation
 *                             values.
 */
Path find_path_multiscale(const Array &z,
                          glm::ivec2   ij_start,
                          glm::ivec2   ij_end,
                          glm::vec4    bbox,
                          int          n_levels = 3,
                          int          corridor_radius = 4,
                          float        corridor_decay = 1.f,
                          float        elevation_ratio = 0.1f,
                          float        distance_exponent = 2.f,
                          float        upward_penalization = 1.f,
                          const Array *p_mask_nogo = nullptr,
                          bool         use_astar = false,
                          bool         smooth_path = false);

} // namespace hmap
