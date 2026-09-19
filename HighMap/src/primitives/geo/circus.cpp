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
#include "highmap/primitives/geo.hpp"

namespace hmap
{

Array circus(glm::ivec2   shape,
             float        radius,
             float        angle,
             float        exit_width,
             float        exit_depth,
             float        center_depth,
             float        ridge_height,
             float        ridge_width,
             float        outer_falloff,
             const Array *p_noise_r,
             glm::vec2    center,
             glm::vec4    bbox,
             Array       *p_mask)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array(shape);

  // --- Setup parameters

  const float rad = std::max(1e-4f, radius);
  const float sigma_w = std::max(1e-4f, ridge_width);
  const float alpha = angle / 180.f * static_cast<float>(M_PI);
  const float exit_w = std::max(1e-3f, exit_width);
  const float falloff_p = std::max(0.1f, outer_falloff);

  Array z = Array(shape);
  Array mask;
  if (p_mask) mask = Array(shape, 0.f);

  // --- Grid coordinates

  bool               endpoint = false;
  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, bbox, endpoint);

  // --- Main generation loop

  for (int j = 0; j < shape.y; j++)
  {
    for (int i = 0; i < shape.x; i++)
    {
      float dx = x[i] - center.x;
      float dy = y[j] - center.y;

      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
      float dist = std::hypot(dx, dy);
      float r = std::max(0.f, dist / rad + dr);

      // polar angle with respect to exit direction
      float theta = std::atan2(dy, dx) - alpha;
      // wrap to [-pi, pi]
      theta = std::atan2(std::sin(theta), std::cos(theta));

      // normalized angle in [0, 1]: 0 at exit (theta = 0), 1 at opposite ridge
      // (theta = +/- pi)
      float norm_angle = 0.5f * (1.f - std::cos(theta));

      // exit opening factor: 1 at exit direction, smoothly falls to 0 away from
      // exit
      float exit_factor = std::exp(-0.5f * (theta * theta) / (exit_w * exit_w));

      // ridge peak elevation along current direction
      // ridge is highest opposite the exit (theta = +/- pi) and dips down to
      // exit_depth at the exit (theta = 0)
      float ridge_rim_z = lerp(exit_depth,
                               ridge_height,
                               norm_angle * (1.f - exit_factor));

      // radial envelope value
      float val = 0.f;

      if (r <= 1.f)
      {
        // inside basin: smoothstep transition from center_depth to ridge_rim_z
        float s = smoothstep3(r);
        val = lerp(center_depth, ridge_rim_z, s);
      }
      else
      {
        // outside basin ridge: smooth falloff towards exit_depth
        float dr_outer = (r - 1.f) / sigma_w;
        float falloff = 1.f / (1.f + std::pow(dr_outer, falloff_p));
        val = lerp(exit_depth, ridge_rim_z, falloff);
      }

      z(i, j) = val;

      if (p_mask)
      {
        // compute smooth mask of the circus shape
        float m = 0.f;
        if (r <= 1.f)
          m = 1.f;
        else
        {
          float dr_outer = (r - 1.f) / sigma_w;
          m = std::max(0.f, 1.f - smoothstep3(std::min(1.f, dr_outer)));
        }
        mask(i, j) = m;
      }
    }
  }

  if (p_mask) *p_mask = mask;

  return z;
}

} // namespace hmap
