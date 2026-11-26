/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/shortest_path.hpp"
#include "highmap/terrain_tri_mesh.hpp"

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
                            bool           favor_sinks,
                            bool           favor_domain_center)
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

  if (favor_domain_center)
  {
    Array mask_nogo(z.shape, 0.f);

    for (int j = 0; j < z.shape.y; ++j)
    {
      for (int i = 0; i < z.shape.x; ++i)
      {
        float di = (z.shape.x > 1)
                       ? std::min(i, z.shape.x - 1 - i) / float(z.shape.x - 1)
                       : 0.f;
        float dj = (z.shape.y > 1)
                       ? std::min(j, z.shape.y - 1 - j) / float(z.shape.y - 1)
                       : 0.f;
        float d_border = std::min(di, dj);
        d_border = std::clamp(d_border, 0.f, 0.5f);

        // normalized penalty: 1 at boundary, 0 at center
        mask_nogo(i, j) = 1.f - 2.f * d_border;
      }
    }

    return find_path_dijkstra(z,
                              start_pt,
                              end_pt,
                              UvRect::unit(),
                              dijk_elevation_ratio,
                              dijk_distance_exponent,
                              dijk_upward_penalization,
                              &mask_nogo);
  }

  return find_path_dijkstra(z,
                            start_pt,
                            end_pt,
                            UvRect::unit(),
                            dijk_elevation_ratio,
                            dijk_distance_exponent,
                            dijk_upward_penalization);
}

Path find_cut_path_dijkstra(const Array   &z,
                            DomainBoundary start,
                            DomainBoundary end,
                            std::uint32_t  seed,
                            size_t         control_points_count,
                            bool           favor_boundary_center,
                            bool           favor_lower_elevation,
                            bool           favor_sinks,
                            bool           favor_domain_center,
                            const UvTile  &region,
                            ICommunicator &comm)
{
  TerrainTriMesh mesh = hmap::generate_terrain_tri_mesh_from_heightmap_random(
      z,
      control_points_count,
      seed,
      region,
      &comm);

  if (comm.rank() != 0 || mesh.size() == 0) return Path();

  size_t start_idx = mesh.pick_mesh_boundary_vertex(start,
                                                    seed,
                                                    favor_boundary_center,
                                                    favor_lower_elevation,
                                                    favor_sinks);

  // use seed + 1 for selecting end vertex to avoid identical index selection on
  // same distributions
  size_t end_idx = mesh.pick_mesh_boundary_vertex(end,
                                                  seed + 1,
                                                  favor_boundary_center,
                                                  favor_lower_elevation,
                                                  favor_sinks);

  float domain_center_weight = favor_domain_center ? 1.f : 0.f;

  std::vector<size_t> vertex_indices = mesh.shortest_path(start_idx,
                                                          end_idx,
                                                          true,
                                                          1.f,
                                                          domain_center_weight);

  if (vertex_indices.empty())
  {
    vertex_indices = {start_idx, end_idx};
  }

  const auto        &pts = mesh.get_points();
  std::vector<Point> path_points;
  path_points.reserve(vertex_indices.size());

  for (size_t v : vertex_indices)
  {
    const auto &p = pts[v];
    path_points.emplace_back(p.x, p.y, p.z);
  }

  return Path(path_points);
}

Path find_cut_path_midpoint(const Array   &z,
                            DomainBoundary start,
                            DomainBoundary end,
                            std::uint32_t  seed,
                            float          offset_ratio,
                            int            steps,
                            bool           favor_boundary_center,
                            bool           favor_lower_elevation,
                            bool           favor_sinks,
                            const UvTile  &region,
                            ICommunicator &comm)
{
  if (!validate_non_empty(z)) return Path();
  if (comm.rank() != 0) return Path();

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
                            region.uv_rect,
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
                              UvRect::unit(),
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
