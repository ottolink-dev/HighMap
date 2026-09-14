/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/shortest_path.hpp"

namespace hmap
{

Path find_cut_path_dijkstra(const Array   &z,
                            DomainBoundary start,
                            DomainBoundary end,
                            float          dijk_elevation_ratio,
                            float          dijk_distance_exponent,
                            float          dijk_upward_penalization,
                            std::uint32_t  seed,
                            bool           favor_boundary_center,
                            bool           favor_lower_elevation,
                            bool           favor_sinks)
{
  if (!validate_non_empty(z)) return Path();

  // --- pick start and end cells

  glm::ivec2 start_pt = pick_boundary_cell(z,
                                           start,
                                           seed,
                                           favor_boundary_center,
                                           favor_lower_elevation,
                                           favor_sinks);
  glm::ivec2 end_pt = pick_boundary_cell(z,
                                         end,
                                         seed,
                                         favor_boundary_center,
                                         favor_lower_elevation,
                                         favor_sinks);

  // --- find cut path

  return find_path_dijkstra(z,
                            start_pt,
                            end_pt,
                            glm::vec4(0.f, 1.f, 0.f, 1.f),
                            dijk_elevation_ratio,
                            dijk_distance_exponent,
                            dijk_upward_penalization);
}

Path find_cut_path_midpoint(const Array   &z,
                            DomainBoundary start,
                            DomainBoundary end,
                            std::uint32_t  seed,
                            float          offset_ratio,
                            int            steps,
                            bool           favor_boundary_center,
                            bool           favor_lower_elevation,
                            bool           favor_sinks)
{
  if (!validate_non_empty(z)) return Path();

  // --- pick start and end cells

  glm::ivec2 start_pt = pick_boundary_cell(z,
                                           start,
                                           seed,
                                           favor_boundary_center,
                                           favor_lower_elevation,
                                           favor_sinks);
  glm::ivec2 end_pt = pick_boundary_cell(z,
                                         end,
                                         seed,
                                         favor_boundary_center,
                                         favor_lower_elevation,
                                         favor_sinks);

  // --- find cut path

  const int max_it = 0; // => autoset by algo

  return find_path_midpoint(z,
                            start_pt,
                            end_pt,
                            glm::vec4(0.f, 1.f, 0.f, 1.f),
                            offset_ratio,
                            max_it,
                            steps);
}

Path find_cut_path_multiscale(const Array   &z,
                              DomainBoundary start,
                              DomainBoundary end,
                              int            n_levels,
                              int            corridor_radius,
                              float          corridor_decay,
                              float          dijk_elevation_ratio,
                              float          dijk_distance_exponent,
                              float          dijk_upward_penalization,
                              std::uint32_t  seed,
                              bool           favor_boundary_center,
                              bool           favor_lower_elevation,
                              bool           favor_sinks,
                              bool           use_astar,
                              bool           smooth_path)
{
  if (!validate_non_empty(z)) return Path();

  // --- pick start and end cells

  glm::ivec2 start_pt = pick_boundary_cell(z,
                                           start,
                                           seed,
                                           favor_boundary_center,
                                           favor_lower_elevation,
                                           favor_sinks);
  glm::ivec2 end_pt = pick_boundary_cell(z,
                                         end,
                                         seed,
                                         favor_boundary_center,
                                         favor_lower_elevation,
                                         favor_sinks);

  // --- find multiscale cut path

  return find_path_multiscale(z,
                              start_pt,
                              end_pt,
                              glm::vec4(0.f, 1.f, 0.f, 1.f),
                              n_levels,
                              corridor_radius,
                              corridor_decay,
                              dijk_elevation_ratio,
                              dijk_distance_exponent,
                              dijk_upward_penalization,
                              nullptr,
                              use_astar,
                              smooth_path);
}

} // namespace hmap
