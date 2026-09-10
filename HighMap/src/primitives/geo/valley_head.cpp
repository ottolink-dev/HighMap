/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"

namespace hmap
{

Array valley_head(glm::ivec2   shape,
                  float        angle,
                  float        max_depth,
                  float        min_width,
                  float        max_width,
                  float        depth_length,
                  float        width_length,
                  float        development_power,
                  float        profile_power,
                  const Array *p_noise_offset,
                  const Array *p_noise_r,
                  glm::vec2    center,
                  glm::vec4    bbox,
                  Array       *p_mask)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_offset && !validate_same_shape(shape, *p_noise_offset))
    return Array(shape);
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array(shape);

  // --- Setup

  const float alpha = angle / 180.f * M_PI;
  const float ca = std::cos(alpha);
  const float sa = std::sin(alpha);

  // clamping parameters for numerical stability
  const float l_d = std::max(1e-6f, depth_length);
  const float l_w = std::max(1e-6f, width_length);
  const float w_min = std::max(1e-6f, min_width);
  const float w_max = std::max(w_min, max_width);
  const float dev_power = std::max(1e-3f, development_power);
  const float prof_power = std::max(1e-3f, profile_power);

  Array z = Array(shape, 0.f);
  Array mask;
  if (p_mask) mask = Array(shape, 0.f);

  // --- Grid coordinates

  bool               endpoint = false;
  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, bbox, endpoint);

  // --- Generate deformation

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      float dx = x[i] - center.x;
      float dy = y[j] - center.y;

      // valley-local coordinates
      // y_local: downstream along valley direction
      // x_local: cross-valley perpendicular
      float x_local = dx * ca - dy * sa;
      float y_local = dx * sa + dy * ca;

      // upstream of valley head: no incision
      if (y_local <= 0.f)
      {
        z(i, j) = 0.f;
        if (p_mask) mask(i, j) = 0.f;
        continue;
      }

      // longitudinal development
      float depth = max_depth *
                    (1.0f - std::exp(-std::pow(y_local / l_d, dev_power)));
      float width = w_min +
                    (w_max - w_min) *
                        (1.0f - std::exp(-std::pow(y_local / l_w, dev_power)));

      float ds = p_noise_offset ? (*p_noise_offset)(i, j) : 0.f;
      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;

      float u = std::abs(x_local / width + ds);
      u = std::max(0.f, u + dr);

      if (u < 1.0f)
      {
        float profile = 1.0f - std::pow(u, prof_power);
        float falloff = 1.f - smoothstep3(u); // lateral
        float factor = profile * falloff;

        z(i, j) = -depth * factor;
        if (p_mask) mask(i, j) = factor;
      }
      else
      {
        z(i, j) = 0.f;
        if (p_mask) mask(i, j) = 0.f;
      }
    }

  if (p_mask) *p_mask = mask;

  return z;
}

} // namespace hmap
