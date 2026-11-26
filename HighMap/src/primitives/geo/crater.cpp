/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/primitives/geo.hpp"

namespace hmap
{

float helper_crater_radial_function(float r,
                                    float inner_depth,
                                    float inner_exp,
                                    float lip_height,
                                    float lip_extent,
                                    float lip_exp,
                                    int   n_terraces,
                                    float terrace_extent,
                                    float terrace_exp,
                                    float terrace_persistence)
{
  // --- Inner crater

  if (r <= 1.f)
  {
    // crater base inner profile
    float v = std::pow(r, inner_exp);
    v = smoothstep3(v);
    v = v * (inner_depth + lip_height) - inner_depth;

    // add inner terraces if any

    // curent terrace index (oriented from lip to crater center)
    float inner = (1.f - r) / terrace_extent;
    int   k = int(inner);

    if (k >= 1 && k < n_terraces + 1)
    {
      // local coordinate inside current terrace
      float u = inner - k;
      float amp = 0.05f * inner_depth * std::pow(terrace_persistence, k - 1);
      v += amp * power_curve(u, terrace_exp, 1.f);
    }

    return v;
  }

  // --- Outside all lips

  float outer = r - 1.f;

  if (outer >= lip_extent) return 0.f;

  // --- Outter lips
  float u = outer / lip_extent;
  float v = smoothstep3(1.f - u);
  return lip_height * std::pow(v, lip_exp);
}

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
                            float terrace_persistence)
{
  float vp = helper_crater_radial_function(r,
                                           inner_depth,
                                           std::max(1.f, ap * inner_exp),
                                           ap * lip_height,
                                           ap * lip_extent,
                                           lip_exp,
                                           n_terraces,
                                           ap * terrace_extent,
                                           terrace_exp,
                                           terrace_persistence);

  float vm = helper_crater_radial_function(r,
                                           inner_depth,
                                           std::max(1.f, am * inner_exp),
                                           am * lip_height,
                                           lip_extent,
                                           lip_exp,
                                           n_terraces,
                                           terrace_extent,
                                           terrace_exp,
                                           terrace_persistence);

  // polar angle w/ respect to main direction (wrapped to [-pi, pi])
  float theta = std::atan2(dy, dx) - alpha;

  // remap angle to [0..1] and use it as a blending factor
  theta = std::atan2(std::sin(theta), std::cos(theta));
  float t = 0.5f * (1.f - std::cos(theta));
  float val = lerp(vp, vm, t);

  // --- Add central peak

  if (central_peak_height != 0.f && r <= central_peak_extent)
  {
    float v = smoothstep3(1.f - r / central_peak_extent);
    val += v * central_peak_height;
  }

  return val;
}

Array crater(glm::ivec2    shape,
             float         radius,
             glm::vec2     center,
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
             Array        *p_inner_crater_mask)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array(shape);

  const float alpha = angle / 180.f * M_PI;
  const float ap = 1.f + asym_ratio;
  const float am = std::max(0.f, 1.f - asym_ratio);

  Array z(shape);
  Array crater_mask(shape, 1.f);
  Array inner_crater_mask(shape);

  // grid
  bool               endpoint = false;
  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, uv_rect, endpoint);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      // --- Compute radius

      float dx = x[i] - center.x;
      float dy = y[j] - center.y;
      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
      float r = dr + std::sqrt(dx * dx + dy * dy) / radius;
      r = std::max(r, 0.f);

      // add sentinel value for the outer region of the crater
      inner_crater_mask(i, j) = r < 1.f ? r : -1.f;

      if (r > 1.f)
      {
        float rm = std::min(1.f, (r - 1.f) / lip_extent);
        crater_mask(i, j) = 1.f - smoothstep3(rm);
      }

      // --- Main profile

      z(i, j) = helper_crater_profile(r,
                                      dx,
                                      dy,
                                      alpha,
                                      ap,
                                      am,
                                      inner_depth,
                                      inner_exp,
                                      lip_height,
                                      lip_extent,
                                      lip_exp,
                                      central_peak_height,
                                      central_peak_extent,
                                      n_terraces,
                                      terrace_extent,
                                      terrace_exp,
                                      terrace_persistence);
    }

  if (p_crater_mask) *p_crater_mask = std::move(crater_mask);
  if (p_inner_crater_mask) *p_inner_crater_mask = std::move(inner_crater_mask);

  return z;
}

} // namespace hmap
