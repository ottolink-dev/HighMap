/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap
{

Array bulkify(const Array         &z,
              const PrimitiveType &primitive_type,
              float                amp,
              const Array         *p_noise_x,
              const Array         *p_noise_y,
              glm::vec2            center,
              const UvRect        &uv_rect)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return Array(z.shape);
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return Array(z.shape);

  return z + amp * get_primitive_base(primitive_type,
                                      z.shape,
                                      p_noise_x,
                                      p_noise_y,
                                      center,
                                      uv_rect);
}

} // namespace hmap
