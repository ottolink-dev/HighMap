/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>

#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

void flatten_around_points_terrain_aware(Array                    &z,
                                         const std::vector<float> &xp,
                                         const std::vector<float> &yp,
                                         const std::vector<float> &zp,
                                         const std::vector<float> &radius,
                                         float                     slope_weight,
                                         float         uphill_weight,
                                         float         max_height_above,
                                         float         terrain_follow,
                                         float         strength,
                                         RadialProfile radial_profile,
                                         float         profile_param,
                                         Array        *p_mask,
                                         glm::vec4     bbox)
{
  if (!validate_non_empty(z)) return;
  if (!validate_non_empty(xp)) return;
  if (!validate_non_empty(yp)) return;
  if (!validate_non_empty(zp)) return;
  if (!validate_non_empty(radius)) return;

  if (xp.size() != yp.size() || xp.size() != zp.size() ||
      xp.size() != radius.size())
    return;

  if (strength <= 0.f) return;

  const int nx = z.shape.x;
  const int ny = z.shape.y;

  // cell size in bbox units (guard against single-cell axes)
  const float dx = nx > 1 ? (bbox.y - bbox.x) / float(nx - 1) : 1.f;
  const float dy = ny > 1 ? (bbox.w - bbox.z) / float(ny - 1) : 1.f;

  auto profile_fct = get_radial_profile_function(radial_profile, profile_param);

  // smoothed terrain (large-scale trend), computed once
  const float max_radius = *std::max_element(radius.begin(), radius.end());
  const int   ir = std::max(
      1,
      int(max_radius / std::min(std::abs(dx), std::abs(dy)) / 2.f));
  Array z_smooth = z;
  smooth_cpulse(z_smooth, ir);

  // per-cell accumulators
  Array weight_sum(z.shape, 0.f);
  Array z_weighted_sum(z.shape, 0.f);
  Array falloff_sum(z.shape, 0.f);

  // cost distance map, reset after each seed (only on visited cells)
  const float inf = std::numeric_limits<float>::max();
  Array       cost(z.shape, inf);

  using Node = std::pair<float, int>; // (cost, linear index)

  // 8-connected neighborhood
  const int di[8] = {1, -1, 0, 0, 1, 1, -1, -1};
  const int dj[8] = {0, 0, 1, -1, 1, -1, 1, -1};

  for (size_t k = 0; k < xp.size(); ++k)
  {
    const float r = radius[k];

    if (r <= 0.f) continue;

    // seed cell
    const int si = int(std::round((xp[k] - bbox.x) / dx));
    const int sj = int(std::round((yp[k] - bbox.z) / dy));

    if (si < 0 || si >= nx || sj < 0 || sj >= ny) continue;

    const float z_seed_smooth = z_smooth(si, sj);
    const float z_limit = max_height_above >= 0.f ? zp[k] + max_height_above
                                                  : inf;

    // --- Dijkstra from the seed

    std::priority_queue<Node, std::vector<Node>, std::greater<Node>> queue;
    std::vector<int>                                                 visited;

    cost(si, sj) = 0.f;
    queue.push({0.f, sj * nx + si});
    visited.push_back(sj * nx + si);

    while (!queue.empty())
    {
      const auto [c, idx] = queue.top();
      queue.pop();

      if (c > cost.vector[idx]) continue;

      const int ci = idx % nx;
      const int cj = idx / nx;

      // --- Accumulate this cell's contribution

      const float u = c / r;
      const float falloff = 1.f - profile_fct(std::clamp(u, 0.f, 1.f));
      const float weight = falloff * r;

      // target: seed elevation + part of the smoothed relief
      const float z_target = zp[k] + terrain_follow *
                                         (z_smooth(ci, cj) - z_seed_smooth);

      weight_sum(ci, cj) += weight;
      z_weighted_sum(ci, cj) += weight * z_target;
      falloff_sum(ci, cj) += falloff;

      // --- Expand to neighbors

      for (int n = 0; n < 8; ++n)
      {
        const int ni = ci + di[n];
        const int nj = cj + dj[n];

        if (ni < 0 || ni >= nx || nj < 0 || nj >= ny) continue;

        const float zn = z(ni, nj);

        if (zn > z_limit) continue;

        // geometric extent limit
        const float ex = bbox.x + float(ni) * dx - xp[k];
        const float ey = bbox.z + float(nj) * dy - yp[k];

        if (ex * ex + ey * ey > r * r) continue;

        const float step_len = std::hypot(float(di[n]) * dx, float(dj[n]) * dy);
        const float slope = std::abs(zn - z(ci, cj)) / step_len;
        const float climb = std::max(0.f, zn - zp[k]) / r;

        const float step_cost = step_len * (1.f + slope_weight * slope +
                                            uphill_weight * climb);
        const float nc = c + step_cost;

        if (nc > r) continue;

        const int nidx = nj * nx + ni;

        if (nc < cost.vector[nidx])
        {
          if (cost.vector[nidx] == inf) visited.push_back(nidx);

          cost.vector[nidx] = nc;
          queue.push({nc, nidx});
        }
      }
    }

    for (int idx : visited)
      cost.vector[idx] = inf;
  }

  // --- Blend toward the weighted target elevation

  if (p_mask) *p_mask = Array(z.shape, 0.f);

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      if (weight_sum(i, j) <= 0.f) continue;

      const float z_ref = z_weighted_sum(i, j) / weight_sum(i, j);
      const float w = std::clamp(strength * falloff_sum(i, j), 0.f, 1.f);

      z(i, j) += w * (z_ref - z(i, j));

      if (p_mask) (*p_mask)(i, j) = w;
    }
}

void flatten_around_points_terrain_aware(Array                    &z,
                                         const std::vector<float> &xp,
                                         const std::vector<float> &yp,
                                         const std::vector<float> &zp,
                                         float                     radius,
                                         float                     slope_weight,
                                         float         uphill_weight,
                                         float         max_height_above,
                                         float         terrain_follow,
                                         float         strength,
                                         RadialProfile radial_profile,
                                         float         profile_param,
                                         Array        *p_mask,
                                         glm::vec4     bbox)
{
  const std::vector<float> radii(xp.size(), radius);

  flatten_around_points_terrain_aware(z,
                                      xp,
                                      yp,
                                      zp,
                                      radii,
                                      slope_weight,
                                      uphill_weight,
                                      max_height_above,
                                      terrain_follow,
                                      strength,
                                      radial_profile,
                                      profile_param,
                                      p_mask,
                                      bbox);
}

} // namespace hmap
