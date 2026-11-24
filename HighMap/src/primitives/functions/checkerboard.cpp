/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdlib>
#include <functional>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap
{

Array checkerboard(glm::ivec2    shape,
                   glm::vec2     kw,
                   const Array  *p_noise_x,
                   const Array  *p_noise_y,
                   const UvRect &uv_rect)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  auto lambda = [&kw](float x, float y, float)
  {
    return std::abs(std::abs((int)std::floor(x) % 2) -
                    std::abs((int)std::floor(y) % 2));
  };

  fill_array_using_xy_function(array,
                               uv_rect,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               lambda);
  return array;
}

} // namespace hmap
