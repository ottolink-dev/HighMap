/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/primitives/geo.hpp"

namespace hmap
{

Array caldera(glm::ivec2    shape,
              float         radius,
              float         sigma_inner,
              float         sigma_outer,
              float         z_bottom,
              const Array  *p_noise,
              float         noise_r_amp,
              float         noise_z_ratio,
              glm::vec2     center,
              const UvRect &uv_rect)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise && !validate_same_shape(shape, *p_noise)) return Array(shape);

  Array z = Array(shape);

  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, uv_rect, false);

  float si2 = sigma_inner * sigma_inner;
  float so2 = sigma_outer * sigma_outer;

  if (p_noise)
  {
    for (int j = 0; j < z.shape.y; j++)
      for (int i = 0; i < z.shape.x; i++)
      {
        float dx = x[i] - center.x;
        float dy = y[j] - center.y;
        float r = std::hypot(dx, dy) - radius;

        r += noise_r_amp * (2 * (*p_noise)(i, j) - 1);

        if (r < 0.f)
          z(i, j) = z_bottom + std::exp(-0.5f * r * r / si2) * (1 - z_bottom);
        else
          z(i, j) = 1.f / (1.f + r * r / so2);

        z(i, j) *= 1.f + noise_z_ratio * (2.f * (*p_noise)(i, j) - 1.f);
      }
  }
  else
  {
    for (int j = 0; j < z.shape.y; j++)
      for (int i = 0; i < z.shape.x; i++)
      {
        float dx = x[i] - center.x;
        float dy = y[j] - center.y;
        float r = std::hypot(dx, dy) - radius;

        if (r < 0.f)
          z(i, j) = z_bottom + std::exp(-0.5f * r * r / si2) * (1 - z_bottom);
        else
          z(i, j) = 1.f / (1.f + r * r / so2);
      }
  }

  return z;
}

Array caldera(glm::ivec2    shape,
              float         radius,
              float         sigma_inner,
              float         sigma_outer,
              float         z_bottom,
              glm::vec2     center,
              const UvRect &uv_rect)
{
  if (!validate_shape(shape)) return Array();

  Array noise = constant(shape, 0.f);
  Array z = caldera(shape,
                    radius,
                    sigma_inner,
                    sigma_outer,
                    z_bottom,
                    nullptr,
                    0.f,
                    0.f,
                    center,
                    uv_rect);
  return z;
}

} // namespace hmap
