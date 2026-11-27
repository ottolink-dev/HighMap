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
#include "highmap/primitives/functions.hpp"

namespace hmap
{

Array polar_shape(glm::ivec2    shape,
                  float         rmin,
                  float         rmax,
                  float         aspect_ratio,
                  float         smoothing_width,
                  bool          square_base,
                  float         angle,
                  float         sector_angle,
                  float         vmin,
                  float         kt_value,
                  float         kr_border,
                  float         kr_border_ratio,
                  const Array  *p_noise_r,
                  const Array  *p_noise_theta,
                  glm::vec2     center,
                  const UvRect &uv_rect)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array(shape);
  if (p_noise_theta && !validate_same_shape(shape, *p_noise_theta))
    return Array(shape);

  const float alpha = angle / 180.f * M_PI;
  const float ca = std::cos(alpha);
  const float sa = std::sin(alpha);
  const float sector_alpha = sector_angle / 180.f * M_PI;
  const float ax = 1.f;
  const float ay = ax / aspect_ratio;

  Array z(shape);
  Array shape_mask(shape);

  // --- Grid

  bool               endpoint = false;
  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, uv_rect, endpoint);

  // --- For each cell loop

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      // --- Compute radius

      float dxc = x[i] - center.x;
      float dyc = y[j] - center.y;

      float dx = (ca * dxc + sa * dyc) / ax;
      float dy = (sa * dxc - ca * dyc) / ay;

      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
      float r = dr + std::sqrt(dx * dx + dy * dy);

      r = std::max(r, 0.f);

      // --- Angle(s)

      // polar angle w/ respect to main direction (in [-pi, pi])
      float theta = std::atan2(dy, dx);

      float dtheta = p_noise_theta ? (*p_noise_theta)(i, j) : 0.f;
      theta += dtheta;

      // remap angle to [0..1] and use it as a blending factor
      theta = std::atan2(std::sin(theta), std::cos(theta));
      float t = 0.5f * (1.f - std::cos(kt_value * theta));

      // --- Value

      float theta_smoothing_width = smoothing_width ? smoothing_width / r : 0.f;

      float value = trapeze_smooth(theta,
                                   -M_PI + 0.5f * sector_alpha,
                                   M_PI - 0.5f * sector_alpha,
                                   theta_smoothing_width);
      value *= lerp(1.f, vmin, t);

      float r1 = rmin;
      float r2 = lerp(rmax,
                      rmin,
                      kr_border_ratio * std::cos(kr_border * theta));

      if (square_base)
      {
        float a = abs_smooth(std::cos(theta) - std::sin(theta), 0.1f);
        float b = abs_smooth(std::cos(theta) + std::sin(theta), 0.1f);
        float denom = a + b;
        r1 /= denom;
        r2 /= denom;
      }

      z(i, j) = value * trapeze_smooth(r, r1, r2, smoothing_width);
    }

  return z;
}

} // namespace hmap
