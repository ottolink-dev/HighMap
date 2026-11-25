/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

void swirl(Array        &dx,
           Array        &dy,
           float         amplitude,
           float         exponent,
           const Array  *p_noise,
           const UvRect &uv_rect)
{
  if (!validate_non_empty(dx) || !validate_same_shape(dx, dy)) return;
  if (p_noise && !validate_same_shape(dx, *p_noise)) return;

  const glm::ivec2 shape = dx.shape;

  glm::vec2 shift = {uv_rect.u_min, uv_rect.v_min};
  glm::vec2 scale = {uv_rect.u_max - uv_rect.u_min,
                     uv_rect.v_max - uv_rect.v_min};

  float xc = (0.5f - shift.x) / scale.x;
  float yc = (0.5f - shift.y) / scale.y;

  if (!p_noise)
  {
    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        float x = (float)i / (float)(shape.x - 1) * scale.x - xc;
        float y = (float)j / (float)(shape.y - 1) * scale.y - yc;
        float re = std::pow(2.f * (x * x + y * y), 0.5f * exponent);
        float theta = 2.f * M_PI * amplitude * re;

        dx(i, j) = std::cos(theta);
        dy(i, j) = std::sin(theta);
      }
  }
  else
  {
    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        float x = (float)i / (float)(shape.x - 1) * scale.x - xc;
        float y = (float)j / (float)(shape.y - 1) * scale.y - yc;
        float re = std::pow(2.f * (x * x + y * y), 0.5f * exponent);
        float theta = 2.f * M_PI * (amplitude * re + (*p_noise)(i, j));

        dx(i, j) = std::cos(theta);
        dy(i, j) = std::sin(theta);
      }
  }
}

} // namespace hmap
