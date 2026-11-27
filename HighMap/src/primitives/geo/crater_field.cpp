/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <random>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/kd_tree.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"
#include "highmap/primitives/geo.hpp"

namespace hmap
{

float helper_crater_profile(float r,
                            float dx,
                            float dy,
                            float alpha,
                            float ap,
                            float am,
                            float inner_depth,
                            float inner_exp,
                            float lip_height,
                            float lip_extent,
                            float lip_exp,
                            float central_peak_height,
                            float central_peak_extent,
                            int   n_terraces,
                            float terrace_extent,
                            float terrace_exp,
                            float terrace_persistence);

Array crater_field(glm::ivec2                    shape,
                   const std::vector<glm::vec2> &centers,
                   const std::vector<float>     &radii,
                   float                         angle,
                   float                         inner_depth,
                   float                         inner_exp,
                   float                         lip_height,
                   float                         lip_extent,
                   float                         lip_exp,
                   float                         asym_ratio,
                   float                         central_peak_height,
                   float                         central_peak_extent,
                   int                           n_terraces,
                   float                         terrace_extent,
                   float                         terrace_exp,
                   float                         terrace_persistence,
                   const Array                  *p_noise_r,
                   const UvRect                 &uv_rect,
                   Array                        *p_crater_mask,
                   Array                        *p_inner_crater_mask,
                   bool                          scale_elevation_with_radius)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array(shape);

  Array z(shape);
  Array crater_mask(shape, 0.f);
  Array inner_crater_mask(shape, -1.f);

  if (centers.size() != radii.size())
    hmap::log::warn("crater_field: centers size ({}) != radii size ({}), "
                    "using min",
                    centers.size(),
                    radii.size());

  const size_t n_craters = std::min(centers.size(), radii.size());

  if (n_craters == 0)
  {
    hmap::log::warn("crater_field: no craters provided, returning empty field");

    if (p_crater_mask) *p_crater_mask = std::move(crater_mask);
    if (p_inner_crater_mask)
      *p_inner_crater_mask = std::move(inner_crater_mask);
    return z;
  }

  // --- Build KD-tree over crater centers

  std::vector<float> x_c(n_craters);
  std::vector<float> y_c(n_craters);
  float              max_radius = 0.f;

  for (size_t k = 0; k < n_craters; k++)
  {
    x_c[k] = centers[k].x;
    y_c[k] = centers[k].y;
    max_radius = std::max(max_radius, radii[k]);
  }

  KDTree tree(x_c, y_c);

  const float alpha = angle / 180.f * M_PI;
  const float ap = 1.f + asym_ratio;
  const float am = std::max(0.f, 1.f - asym_ratio);
  const float max_lip_scale = 1.f + std::max(ap, 1.f) * lip_extent;
  const float max_influence_radius = max_radius * max_lip_scale;

  // --- Precompute sequential local reference elevation at each crater center

  std::vector<float> z_ref(n_craters, 0.f);

  for (size_t k = 0; k < n_craters; k++)
  {
    for (size_t m = 0; m < k; m++)
    {
      const float radius_m = radii[m];
      if (radius_m <= 0.f) continue;

      const float dx = centers[k].x - centers[m].x;
      const float dy = centers[k].y - centers[m].y;
      const float r = std::sqrt(dx * dx + dy * dy) / radius_m;

      if (r >= max_lip_scale) continue;

      float local_mask = 1.f;
      if (r >= 1.f)
      {
        float rm = std::min(1.f, (r - 1.f) / lip_extent);
        local_mask = 1.f - smoothstep3(rm);
      }

      const float elev_scale = scale_elevation_with_radius ? radius_m : 1.f;

      const float val = helper_crater_profile(r,
                                              dx,
                                              dy,
                                              alpha,
                                              ap,
                                              am,
                                              inner_depth * elev_scale,
                                              inner_exp,
                                              lip_height * elev_scale,
                                              lip_extent,
                                              lip_exp,
                                              central_peak_height * elev_scale,
                                              central_peak_extent,
                                              n_terraces,
                                              terrace_extent,
                                              terrace_exp,
                                              terrace_persistence);

      z_ref[k] = lerp(z_ref[k], z_ref[m], local_mask);
      z_ref[k] += local_mask * val;
    }
  }

  // --- Grid

  bool               endpoint = false;
  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, uv_rect, endpoint);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      auto candidates = tree.radius_search(x[i], y[j], max_influence_radius);

      if (candidates.empty()) continue;

      // sort by crater index to apply impacts sequentially
      std::sort(candidates.begin(),
                candidates.end(),
                [](const auto &a, const auto &b) { return a.first < b.first; });

      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;

      for (const auto &[k, dist_sq] : candidates)
      {
        const float radius = radii[k];
        if (radius <= 0.f) continue;

        float dx = x[i] - centers[k].x;
        float dy = y[j] - centers[k].y;
        float r = dr + std::sqrt(dist_sq) / radius;
        r = std::max(r, 0.f);

        if (r >= max_lip_scale) continue;

        // update masks
        float local_mask = 1.f;

        if (r < 1.f)
        {
          crater_mask(i, j) = 1.f;

          if (inner_crater_mask(i, j) < 0.f || r < inner_crater_mask(i, j))
            inner_crater_mask(i, j) = r;
        }
        else
        {
          float rm = std::min(1.f, (r - 1.f) / lip_extent);
          local_mask = 1.f - smoothstep3(rm);
          crater_mask(i, j) = std::max(crater_mask(i, j), local_mask);
        }

        const float elev_scale = scale_elevation_with_radius ? radius : 1.f;

        float val = helper_crater_profile(r,
                                          dx,
                                          dy,
                                          alpha,
                                          ap,
                                          am,
                                          inner_depth * elev_scale,
                                          inner_exp,
                                          lip_height * elev_scale,
                                          lip_extent,
                                          lip_exp,
                                          central_peak_height * elev_scale,
                                          central_peak_extent,
                                          n_terraces,
                                          terrace_extent,
                                          terrace_exp,
                                          terrace_persistence);

        z(i, j) = lerp(z(i, j), z_ref[k], local_mask);
        z(i, j) += local_mask * val;
      }
    }

  if (p_crater_mask) *p_crater_mask = std::move(crater_mask);
  if (p_inner_crater_mask) *p_inner_crater_mask = std::move(inner_crater_mask);

  return z;
}

Array crater_field(glm::ivec2    shape,
                   uint          seed,
                   int           n_craters,
                   float         r_min,
                   float         r_max,
                   float         angle,
                   float         inner_depth,
                   float         inner_exp,
                   float         lip_height,
                   float         lip_extent,
                   float         lip_exp,
                   float         asym_ratio,
                   float         central_peak_height,
                   float         central_peak_extent,
                   int           n_terraces,
                   float         terrace_extent,
                   float         terrace_exp,
                   float         terrace_persistence,
                   const Array  *p_noise_r,
                   const UvRect &uv_rect,
                   Array        *p_crater_mask,
                   Array        *p_inner_crater_mask,
                   bool          scale_elevation_with_radius)
{
  if (n_craters <= 0)
  {
    hmap::log::warn("crater_field: n_craters ({}) <= 0, returning empty field",
                    n_craters);

    return crater_field(shape,
                        std::vector<glm::vec2>{},
                        std::vector<float>{},
                        angle,
                        inner_depth,
                        inner_exp,
                        lip_height,
                        lip_extent,
                        lip_exp,
                        asym_ratio,
                        central_peak_height,
                        central_peak_extent,
                        n_terraces,
                        terrace_extent,
                        terrace_exp,
                        terrace_persistence,
                        p_noise_r,
                        uv_rect,
                        p_crater_mask,
                        p_inner_crater_mask);
  }

  // sample centers in the unit square [0, 1]^2 so the algorithm is tileable
  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dist_xy(0.f, 1.f);
  std::uniform_real_distribution<float> dist_r(std::min(r_min, r_max),
                                               std::max(r_min, r_max));

  std::vector<glm::vec2> centers(n_craters);
  std::vector<float>     radii(n_craters);

  for (int k = 0; k < n_craters; k++)
  {
    centers[k] = {dist_xy(gen), dist_xy(gen)};
    radii[k] = dist_r(gen);
  }

  return crater_field(shape,
                      centers,
                      radii,
                      angle,
                      inner_depth,
                      inner_exp,
                      lip_height,
                      lip_extent,
                      lip_exp,
                      asym_ratio,
                      central_peak_height,
                      central_peak_extent,
                      n_terraces,
                      terrace_extent,
                      terrace_exp,
                      terrace_persistence,
                      p_noise_r,
                      uv_rect,
                      p_crater_mask,
                      p_inner_crater_mask,
                      scale_elevation_with_radius);
}

} // namespace hmap
