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

  const int nx = z.shape.x;
  const int ny = z.shape.y;

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

  std::vector<int> i_path, j_path;

  find_path_dijkstra(z,
                     glm::ivec2(start_pt.x, start_pt.y),
                     glm::ivec2(end_pt.x, end_pt.y),
                     i_path,
                     j_path,
                     dijk_elevation_ratio,
                     dijk_distance_exponent,
                     dijk_upward_penalization);

  // --- build the output path

  std::vector<float> x, y, v;
  x.reserve(i_path.size());
  y.reserve(i_path.size());
  v.reserve(i_path.size());

  for (size_t k = 0; k < i_path.size(); ++k)
  {
    x.push_back(float(i_path[k]) / float(nx - 1));
    y.push_back(float(j_path[k]) / float(ny - 1));
    v.push_back(z(i_path[k], j_path[k]));
  }

  return Path(x, y, v);
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

  const int nx = z.shape.x;
  const int ny = z.shape.y;

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

  int max_it = 0; // => autoset by algo

  std::vector<glm::ivec2> indices = find_path_midpoint(z,
                                                       start_pt,
                                                       end_pt,
                                                       offset_ratio,
                                                       max_it,
                                                       steps);

  const size_t       npts = indices.size();
  std::vector<float> x, y, v;

  x.reserve(npts);
  y.reserve(npts);
  v.reserve(npts);

  for (const auto &p : indices)
  {
    x.push_back(float(p.x) / float(nx - 1));
    y.push_back(float(p.y) / float(ny - 1));
    v.push_back(z(p));
  }

  auto path = Path(x, y, v);
  return path;
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
