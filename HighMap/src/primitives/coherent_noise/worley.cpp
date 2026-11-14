/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/coherent_noise.hpp"

namespace hmap
{

Array worley_double(glm::ivec2    shape,
                    glm::vec2     kw,
                    std::uint32_t seed,
                    float         ratio,
                    float         k,
                    const Array  *p_ctrl_param,
                    const Array  *p_noise_x,
                    const Array  *p_noise_y,
                    glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  hmap::Array                array = hmap::Array(shape);
  hmap::WorleyDoubleFunction f = hmap::WorleyDoubleFunction(kw, seed, ratio, k);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f);
  return array;
}

} // namespace hmap
