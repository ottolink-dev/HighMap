/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/shortest_path.hpp"

namespace hmap
{

namespace
{

struct SearchNode
{
  float priority;
  float g_cost;
  int   i;
  int   j;

  bool operator>(const SearchNode &other) const
  {
    return priority > other.priority;
  }
};

// --- Helper Functions

std::vector<glm::ivec2> helper_solve_graph_search(
    const Array        &z,
    glm::ivec2          ij_start,
    glm::ivec2          ij_end,
    const Mat<uint8_t> *p_corridor,
    float               elevation_ratio,
    float               distance_exponent,
    float               upward_penalization,
    const Array        *p_mask_nogo,
    bool                use_astar)
{
  glm::ivec2 shape = z.shape;

  // 8-neighbor directions and weights
  const std::vector<int>   di = {-1, 0, 0, 1, -1, -1, 1, 1};
  const std::vector<int>   dj = {0, 1, -1, 0, -1, 1, -1, 1};
  const std::vector<float> cd = {1.f,
                                 1.f,
                                 1.f,
                                 1.f,
                                 float(M_SQRT2),
                                 float(M_SQRT2),
                                 float(M_SQRT2),
                                 float(M_SQRT2)};
  const size_t             nb = di.size();

  Mat<float> g_score(shape, std::numeric_limits<float>::infinity());
  Mat<int>   parent_i(shape, -1);
  Mat<int>   parent_j(shape, -1);

  std::priority_queue<SearchNode,
                      std::vector<SearchNode>,
                      std::greater<SearchNode>>
      pq;

  auto heuristic = [&](int i, int j) -> float
  {
    if (!use_astar) return 0.f;

    // admissible elevation heuristic: minimum cost to reach destination
    // elevation
    float dz = z(ij_end.x, ij_end.y) - z(i, j);
    return elevation_ratio * std::max(0.f, dz);
  };

  g_score(ij_start.x, ij_start.y) = 0.f;
  pq.push({heuristic(ij_start.x, ij_start.y), 0.f, ij_start.x, ij_start.y});

  bool reached_target = false;

  while (!pq.empty())
  {
    SearchNode curr = pq.top();
    pq.pop();

    int i = curr.i;
    int j = curr.j;

    if (curr.g_cost > g_score(i, j)) continue;

    if (i == ij_end.x && j == ij_end.y)
    {
      reached_target = true;
      break;
    }

    for (size_t k = 0; k < nb; k++)
    {
      int p = i + di[k];
      int q = j + dj[k];

      if (p >= 0 && p < shape.x && q >= 0 && q < shape.y)
      {
        // check corridor filter if active
        if (p_corridor && !(*p_corridor)(p, q)) continue;

        // elevation difference contribution
        float dz = (z(i, j) - z(p, q)) * cd[k];
        if (dz < 0.f) dz *= upward_penalization;
        dz = std::abs(dz);

        float step_cost = (1.f - elevation_ratio) *
                          std::pow(dz, distance_exponent);

        // absolute elevation contribution
        step_cost += elevation_ratio *
                     std::max(0.f, cd[k] * (z(p, q) - z(i, j)));

        if (p_mask_nogo) step_cost += 1e5f * (*p_mask_nogo)(p, q);

        float new_g = curr.g_cost + step_cost;

        if (new_g < g_score(p, q))
        {
          g_score(p, q) = new_g;
          parent_i(p, q) = i;
          parent_j(p, q) = j;
          pq.push({new_g + heuristic(p, q), new_g, p, q});
        }
      }
    }
  }

  if (!reached_target) return {};

  // --- Build path backwards

  std::vector<glm::ivec2> path;
  int                     ic = ij_end.x;
  int                     jc = ij_end.y;

  while (ic != ij_start.x || jc != ij_start.y)
  {
    path.push_back({ic, jc});
    int pi = parent_i(ic, jc);
    int pj = parent_j(ic, jc);
    if (pi == -1 || pj == -1) break;
    ic = pi;
    jc = pj;
  }

  path.push_back({ij_start.x, ij_start.y});
  std::reverse(path.begin(), path.end());

  return path;
}

std::vector<glm::ivec2> helper_smooth_grid_path(
    const std::vector<glm::ivec2> &indices,
    glm::ivec2                     shape)
{
  if (indices.size() < 3) return indices;

  // build path with normalized coordinates
  std::vector<float> x, y;
  x.reserve(indices.size());
  y.reserve(indices.size());

  for (const auto &p : indices)
  {
    x.push_back(float(p.x) / float(shape.x - 1));
    y.push_back(float(p.y) / float(shape.y - 1));
  }

  Path path(x, y);
  Path smoothed = smooth(path, 2, 1.f, 0.f);

  // convert smoothed path back to grid indices
  std::vector<glm::ivec2> result;
  result.reserve(smoothed.size());
  result.push_back(indices.front());

  for (size_t k = 1; k < smoothed.size(); ++k)
  {
    Point p0 = smoothed.points[k - 1];
    Point p1 = smoothed.points[k];

    float u0 = p0.x * float(shape.x - 1);
    float v0 = p0.y * float(shape.y - 1);
    float u1 = p1.x * float(shape.x - 1);
    float v1 = p1.y * float(shape.y - 1);

    float seg_len = std::hypot(u1 - u0, v1 - v0);
    int   n_steps = std::max(1, int(std::ceil(seg_len * 2.f)));

    for (int s = 1; s <= n_steps; ++s)
    {
      float t = float(s) / float(n_steps);
      int   sx = int(std::round(u0 + t * (u1 - u0)));
      int   sy = int(std::round(v0 + t * (v1 - v0)));
      sx = std::clamp(sx, 0, shape.x - 1);
      sy = std::clamp(sy, 0, shape.y - 1);

      glm::ivec2 ij(sx, sy);
      if (ij != result.back()) result.push_back(ij);
    }
  }

  if (result.back() != indices.back()) result.push_back(indices.back());

  return result;
}

} // namespace

// --- Main Multiscale Functions

void find_path_multiscale(const Array      &z,
                          glm::ivec2        ij_start,
                          glm::ivec2        ij_end,
                          std::vector<int> &i_path,
                          std::vector<int> &j_path,
                          int               n_levels,
                          int               corridor_radius,
                          float             corridor_decay,
                          float             elevation_ratio,
                          float             distance_exponent,
                          float             upward_penalization,
                          const Array      *p_mask_nogo,
                          bool              use_astar,
                          bool              smooth_path)
{
  i_path.clear();
  j_path.clear();

  std::vector<glm::ivec2> path = find_path_multiscale(z,
                                                      ij_start,
                                                      ij_end,
                                                      n_levels,
                                                      corridor_radius,
                                                      corridor_decay,
                                                      elevation_ratio,
                                                      distance_exponent,
                                                      upward_penalization,
                                                      p_mask_nogo,
                                                      use_astar,
                                                      smooth_path);

  i_path.reserve(path.size());
  j_path.reserve(path.size());

  for (const auto &pt : path)
  {
    i_path.push_back(pt.x);
    j_path.push_back(pt.y);
  }
}

std::vector<glm::ivec2> find_path_multiscale(const Array &z,
                                             glm::ivec2   ij_start,
                                             glm::ivec2   ij_end,
                                             int          n_levels,
                                             int          corridor_radius,
                                             float        corridor_decay,
                                             float        elevation_ratio,
                                             float        distance_exponent,
                                             float        upward_penalization,
                                             const Array *p_mask_nogo,
                                             bool         use_astar,
                                             bool         smooth_path)
{
  if (!validate_non_empty(z)) return {};
  if (p_mask_nogo && !validate_same_shape(z, *p_mask_nogo)) return {};

  const glm::ivec2 &shape_orig = z.shape;

  ij_start.x = std::clamp(ij_start.x, 0, shape_orig.x - 1);
  ij_start.y = std::clamp(ij_start.y, 0, shape_orig.y - 1);
  ij_end.x = std::clamp(ij_end.x, 0, shape_orig.x - 1);
  ij_end.y = std::clamp(ij_end.y, 0, shape_orig.y - 1);

  if (ij_start == ij_end) return {ij_start};

  // --- Determine pyramid level shapes

  n_levels = std::max(1, n_levels);
  std::vector<glm::ivec2> level_shapes;
  level_shapes.push_back(shape_orig);

  for (int l = 1; l < n_levels; ++l)
  {
    int        factor = 1 << l;
    glm::ivec2 sh = {std::max(4, shape_orig.x / factor),
                     std::max(4, shape_orig.y / factor)};
    if (sh == level_shapes.back()) break;
    level_shapes.push_back(sh);
  }

  int actual_levels = int(level_shapes.size());

  // --- Build downsampled pyramid

  std::vector<Array> z_pyramid(actual_levels);
  std::vector<Array> mask_pyramid(actual_levels);

  z_pyramid[0] = z;
  if (p_mask_nogo) mask_pyramid[0] = *p_mask_nogo;

  for (int l = 1; l < actual_levels; ++l)
  {
    z_pyramid[l] = z.resample_to_shape(level_shapes[l]);
    if (p_mask_nogo)
      mask_pyramid[l] = p_mask_nogo->resample_to_shape(level_shapes[l]);
  }

  // --- Coarsest level solve

  int        coarsest = actual_levels - 1;
  glm::ivec2 shape_c = level_shapes[coarsest];

  glm::ivec2 start_c = {
      int(std::round(float(ij_start.x) / float(shape_orig.x - 1) *
                     (shape_c.x - 1))),
      int(std::round(float(ij_start.y) / float(shape_orig.y - 1) *
                     (shape_c.y - 1)))};
  glm::ivec2 end_c = {int(std::round(float(ij_end.x) / float(shape_orig.x - 1) *
                                     (shape_c.x - 1))),
                      int(std::round(float(ij_end.y) / float(shape_orig.y - 1) *
                                     (shape_c.y - 1)))};

  start_c.x = std::clamp(start_c.x, 0, shape_c.x - 1);
  start_c.y = std::clamp(start_c.y, 0, shape_c.y - 1);
  end_c.x = std::clamp(end_c.x, 0, shape_c.x - 1);
  end_c.y = std::clamp(end_c.y, 0, shape_c.y - 1);

  const Array *p_mask_c = p_mask_nogo ? &mask_pyramid[coarsest] : nullptr;

  std::vector<glm::ivec2> current_path = helper_solve_graph_search(
      z_pyramid[coarsest],
      start_c,
      end_c,
      nullptr,
      elevation_ratio,
      distance_exponent,
      upward_penalization,
      p_mask_c,
      use_astar);

  if (current_path.empty())
  {
    // fallback to direct line if obstructed
    current_path = {start_c, end_c};
  }

  // --- Propagate coarse-to-fine

  for (int l = coarsest - 1; l >= 0; --l)
  {
    glm::ivec2 prev_shape = level_shapes[l + 1];
    glm::ivec2 cur_shape = level_shapes[l];

    glm::ivec2 cur_start = (l == 0)
                               ? ij_start
                               : glm::ivec2(
                                     int(std::round(float(ij_start.x) /
                                                    float(shape_orig.x - 1) *
                                                    (cur_shape.x - 1))),
                                     int(std::round(float(ij_start.y) /
                                                    float(shape_orig.y - 1) *
                                                    (cur_shape.y - 1))));

    glm::ivec2 cur_end = (l == 0) ? ij_end
                                  : glm::ivec2(
                                        int(std::round(float(ij_end.x) /
                                                       float(shape_orig.x - 1) *
                                                       (cur_shape.x - 1))),
                                        int(std::round(float(ij_end.y) /
                                                       float(shape_orig.y - 1) *
                                                       (cur_shape.y - 1))));

    cur_start.x = std::clamp(cur_start.x, 0, cur_shape.x - 1);
    cur_start.y = std::clamp(cur_start.y, 0, cur_shape.y - 1);
    cur_end.x = std::clamp(cur_end.x, 0, cur_shape.x - 1);
    cur_end.y = std::clamp(cur_end.y, 0, cur_shape.y - 1);

    // compute corridor radius for this level
    int cur_radius = std::max(
        1,
        int(std::round(float(corridor_radius) *
                       std::pow(corridor_decay, float(l)))));

    // build corridor mask at current resolution
    Mat<uint8_t> corridor(cur_shape, uint8_t(0));

    auto mark_corridor_cell = [&](int cx, int cy)
    {
      for (int dx = -cur_radius; dx <= cur_radius; ++dx)
      {
        for (int dy = -cur_radius; dy <= cur_radius; ++dy)
        {
          if (dx * dx + dy * dy <= cur_radius * cur_radius)
          {
            int px = cx + dx;
            int py = cy + dy;
            if (px >= 0 && px < cur_shape.x && py >= 0 && py < cur_shape.y)
            {
              corridor(px, py) = 1;
            }
          }
        }
      }
    };

    // sample upscaled path densely into corridor
    for (size_t k = 0; k < current_path.size(); ++k)
    {
      float u0 = float(current_path[k].x) / float(prev_shape.x - 1) *
                 (cur_shape.x - 1);
      float v0 = float(current_path[k].y) / float(prev_shape.y - 1) *
                 (cur_shape.y - 1);

      if (k + 1 < current_path.size())
      {
        float u1 = float(current_path[k + 1].x) / float(prev_shape.x - 1) *
                   (cur_shape.x - 1);
        float v1 = float(current_path[k + 1].y) / float(prev_shape.y - 1) *
                   (cur_shape.y - 1);

        float seg_len = std::hypot(u1 - u0, v1 - v0);
        int   steps = std::max(1, int(std::ceil(seg_len * 2.f)));

        for (int s = 0; s <= steps; ++s)
        {
          float t = float(s) / float(steps);
          int   sx = int(std::round(u0 + t * (u1 - u0)));
          int   sy = int(std::round(v0 + t * (v1 - v0)));
          mark_corridor_cell(sx, sy);
        }
      }
      else
      {
        mark_corridor_cell(int(std::round(u0)), int(std::round(v0)));
      }
    }

    mark_corridor_cell(cur_start.x, cur_start.y);
    mark_corridor_cell(cur_end.x, cur_end.y);

    const Array *p_mask_cur = p_mask_nogo ? &mask_pyramid[l] : nullptr;

    std::vector<glm::ivec2> next_path = helper_solve_graph_search(
        z_pyramid[l],
        cur_start,
        cur_end,
        &corridor,
        elevation_ratio,
        distance_exponent,
        upward_penalization,
        p_mask_cur,
        use_astar);

    if (next_path.empty())
    {
      // fallback to unconstrained search on failure
      next_path = helper_solve_graph_search(z_pyramid[l],
                                            cur_start,
                                            cur_end,
                                            nullptr,
                                            elevation_ratio,
                                            distance_exponent,
                                            upward_penalization,
                                            p_mask_cur,
                                            use_astar);
    }

    if (!next_path.empty()) current_path = std::move(next_path);
  }

  // --- Path smoothing

  if (smooth_path && current_path.size() > 2)
  {
    current_path = helper_smooth_grid_path(current_path, shape_orig);
  }

  return current_path;
}

Path find_path_multiscale(const Array &z,
                          glm::ivec2   ij_start,
                          glm::ivec2   ij_end,
                          glm::vec4    bbox,
                          int          n_levels,
                          int          corridor_radius,
                          float        corridor_decay,
                          float        elevation_ratio,
                          float        distance_exponent,
                          float        upward_penalization,
                          const Array *p_mask_nogo,
                          bool         use_astar,
                          bool         smooth_path)
{
  std::vector<glm::ivec2> indices = find_path_multiscale(z,
                                                         ij_start,
                                                         ij_end,
                                                         n_levels,
                                                         corridor_radius,
                                                         corridor_decay,
                                                         elevation_ratio,
                                                         distance_exponent,
                                                         upward_penalization,
                                                         p_mask_nogo,
                                                         use_astar,
                                                         false);

  if (indices.empty()) return Path();

  const int nx = z.shape.x;
  const int ny = z.shape.y;

  Path path;
  path.points.reserve(indices.size());

  for (const auto &p : indices)
  {
    float x = (nx > 1) ? float(p.x) / float(nx - 1) : 0.f;
    float y = (ny > 1) ? float(p.y) / float(ny - 1) : 0.f;

    x = bbox.x + x * (bbox.y - bbox.x);
    y = bbox.z + y * (bbox.w - bbox.z);
    float v = z(p.x, p.y);

    path.add_point(Point(x, y, v));
  }

  if (smooth_path && path.size() > 2)
  {
    path = smooth(path);
  }

  return path;
}

} // namespace hmap
