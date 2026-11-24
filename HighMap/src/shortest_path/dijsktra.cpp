/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/shortest_path.hpp"
#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array.hpp"

#include <format>

namespace hmap
{

void find_path_dijkstra(const Array                   &z,
                        glm::ivec2                     ij_start,
                        std::vector<glm::ivec2>        ij_end_list,
                        std::vector<std::vector<int>> &i_path_list,
                        std::vector<std::vector<int>> &j_path_list,
                        float                          elevation_ratio,
                        float                          distance_exponent,
                        float                          upward_penalization,
                        const Array                   *p_mask_nogo,
                        const Array                   *p_node_weights)
{
  if (!validate_non_empty(z)) return;
  if (p_mask_nogo && !validate_same_shape(z, *p_mask_nogo)) return;
  if (p_node_weights && !validate_same_shape(z, *p_node_weights)) return;

  // https://math.stackexchange.com/questions/3088292

  glm::ivec2 shape = z.shape;

  // neighbors pattern
  const auto  &neighbors = neighborhood::MOORE_8;
  const size_t nb = neighbors.size();

  // working arrays
  Mat<float> distance(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      distance(i, j) = std::numeric_limits<float>::infinity();

  Mat<int> next_idx_i(shape);
  Mat<int> next_idx_j(shape);

  // min-priority queue: stores (distance, packed_cell_index)
  using QueueElement = std::pair<float, int>;
  std::priority_queue<QueueElement,
                      std::vector<QueueElement>,
                      std::greater<QueueElement>>
      queue;

  // --- Dijkstra's algorithm

  int start_idx = ij_start.y * shape.x + ij_start.x;
  distance(ij_start.x, ij_start.y) = 0.f;
  queue.emplace(0.f, start_idx);

  while (!queue.empty())
  {
    auto [curr_dist, idx] = queue.top();
    queue.pop();

    int i = idx % shape.x;
    int j = idx / shape.x;

    if (curr_dist > distance(i, j)) continue;

    // loop over neighbors
    for (size_t k = 0; k < nb; k++)
    {
      int p = i + neighbors[k].offset.x;
      int q = j + neighbors[k].offset.y;

      if ((p >= 0) and (p < shape.x) and (q >= 0) and (q < shape.y))
      {
        // elevation difference contribution (weighted for diagonal
        // directions to avoid artifacts)
        float dz = (z(i, j) - z(p, q)) * neighbors[k].distance;
        if (dz < 0.f) dz *= upward_penalization;
        dz = std::abs(dz);

        float edge_cost = (1.f - elevation_ratio) *
                          std::pow(dz, distance_exponent);

        // absolute elevation contribution (puts the emphasis on
        // going downslope rather than upslope)
        edge_cost += elevation_ratio *
                     std::max(0.f, neighbors[k].distance * (z(p, q) - z(i, j)));

        if (p_mask_nogo) edge_cost += 1e5f * (*p_mask_nogo)(p, q);

        if (p_node_weights) edge_cost *= (*p_node_weights)(p, q);

        float new_dist = curr_dist + edge_cost;

        if (new_dist < distance(p, q))
        {
          distance(p, q) = new_dist;

          next_idx_i(p, q) = i;
          next_idx_j(p, q) = j;

          queue.emplace(new_dist, q * shape.x + p);
        }
      }
    }
  }

  // --- Build path backwards

  i_path_list.clear();
  j_path_list.clear();

  for (auto ij_end : ij_end_list)
  {
    int ic = ij_end.x;
    int jc = ij_end.y;

    std::vector<int> i_path, j_path;

    if (distance(ic, jc) < std::numeric_limits<float>::infinity())
    {
      while ((ic != ij_start.x) or (jc != ij_start.y))
      {
        i_path.push_back(ic);
        j_path.push_back(jc);

        int prev_i = next_idx_i(ic, jc);
        int prev_j = next_idx_j(ic, jc);
        ic = prev_i;
        jc = prev_j;
      }

      i_path.push_back(ic);
      j_path.push_back(jc);

      std::reverse(i_path.begin(), i_path.end());
      std::reverse(j_path.begin(), j_path.end());
    }

    i_path_list.push_back(i_path);
    j_path_list.push_back(j_path);
  }
}

void find_path_dijkstra(const Array      &z,
                        glm::ivec2        ij_start,
                        glm::ivec2        ij_end,
                        std::vector<int> &i_path,
                        std::vector<int> &j_path,
                        float             elevation_ratio,
                        float             distance_exponent,
                        float             upward_penalization,
                        const Array      *p_mask_nogo,
                        const Array      *p_node_weights)
{
  std::vector<std::vector<int>> i_path_list = {i_path};
  std::vector<std::vector<int>> j_path_list = {j_path};

  find_path_dijkstra(z,
                     ij_start,
                     {ij_end},
                     i_path_list,
                     j_path_list,
                     elevation_ratio,
                     distance_exponent,
                     upward_penalization,
                     p_mask_nogo,
                     p_node_weights);

  i_path = i_path_list[0];
  j_path = j_path_list[0];
}

std::vector<glm::ivec2> find_path_dijkstra(const Array &z,
                                           glm::ivec2   ij_start,
                                           glm::ivec2   ij_end,
                                           float        elevation_ratio,
                                           float        distance_exponent,
                                           float        upward_penalization,
                                           const Array *p_mask_nogo,
                                           const Array *p_node_weights)
{
  std::vector<int> ip, jp;
  find_path_dijkstra(z,
                     ij_start,
                     ij_end,
                     ip,
                     jp,
                     elevation_ratio,
                     distance_exponent,
                     upward_penalization,
                     p_mask_nogo,
                     p_node_weights);

  std::vector<glm::ivec2> path;
  path.reserve(ip.size());
  for (size_t k = 0; k < ip.size(); ++k)
    path.emplace_back(ip[k], jp[k]);

  return path;
}

Path find_path_dijkstra(const Array  &z,
                        glm::ivec2    ij_start,
                        glm::ivec2    ij_end,
                        const UvRect &uv_rect,
                        float         elevation_ratio,
                        float         distance_exponent,
                        float         upward_penalization,
                        const Array  *p_mask_nogo,
                        const Array  *p_node_weights)
{
  if (!validate_non_empty(z)) return Path();

  std::vector<int> ip, jp;
  find_path_dijkstra(z,
                     ij_start,
                     ij_end,
                     ip,
                     jp,
                     elevation_ratio,
                     distance_exponent,
                     upward_penalization,
                     p_mask_nogo,
                     p_node_weights);

  std::vector<Point> points;
  points.reserve(ip.size());

  const float lx = uv_rect.u_max - uv_rect.u_min;
  const float ly = uv_rect.v_max - uv_rect.v_min;
  const float denom_x = (z.shape.x > 1) ? float(z.shape.x - 1) : 1.f;
  const float denom_y = (z.shape.y > 1) ? float(z.shape.y - 1) : 1.f;

  for (size_t k = 0; k < ip.size(); ++k)
  {
    float px = (float(ip[k]) / denom_x) * lx + uv_rect.u_min;
    float py = (float(jp[k]) / denom_y) * ly + uv_rect.v_min;
    float pv = z(ip[k], jp[k]);

    points.emplace_back(px, py, pv);
  }

  return Path(points);
}

} // namespace hmap

namespace hmap::va
{

Path find_path_dijkstra(const ComputeMode  &cm,
                        const VirtualArray &z,
                        glm::vec2           start,
                        glm::vec2           end,
                        std::uint32_t       seed,
                        size_t              control_points_count,
                        bool                favor_domain_center)
{
  TerrainTriMesh mesh = va::generate_terrain_tri_mesh_from_heightmap_random(
      cm,
      z,
      control_points_count,
      seed);

  if (mesh.size() == 0) return Path();

  const auto &pts = mesh.get_points();

  // find nearest vertex to start position
  size_t start_idx = 0;
  float  min_dist2_start = std::numeric_limits<float>::max();

  for (size_t i = 0; i < pts.size(); ++i)
  {
    float dx = pts[i].x - start.x;
    float dy = pts[i].y - start.y;
    float d2 = dx * dx + dy * dy;

    if (d2 < min_dist2_start)
    {
      min_dist2_start = d2;
      start_idx = i;
    }
  }

  // find nearest vertex to end position
  size_t end_idx = 0;
  float  min_dist2_end = std::numeric_limits<float>::max();

  for (size_t i = 0; i < pts.size(); ++i)
  {
    float dx = pts[i].x - end.x;
    float dy = pts[i].y - end.y;
    float d2 = dx * dx + dy * dy;

    if (d2 < min_dist2_end)
    {
      min_dist2_end = d2;
      end_idx = i;
    }
  }

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

  std::vector<Point> path_points;
  path_points.reserve(vertex_indices.size());

  for (size_t v : vertex_indices)
  {
    const auto &p = pts[v];
    path_points.emplace_back(p.x, p.y, p.z);
  }

  return Path(path_points);
}

} // namespace hmap::va
