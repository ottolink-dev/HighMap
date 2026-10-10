/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap::gpu
{

Array advection_warp(const Array &z,
                     const Array &advected_field,
                     const Array &dx,
                     const Array &dy,
                     float        advection_length,
                     float        value_persistence,
                     const Array *p_mask)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, advected_field) ||
      !validate_same_shape(z, dx) || !validate_same_shape(z, dy))
    return Array();
  if (p_mask && !validate_same_shape(z, *p_mask)) return Array();

  auto run = clwrapper::Run("advection_warp");

  glm::ivec2 shape = z.shape;
  Array      out = Array(shape);
  Array      mask;

  run.bind_imagef("z", z.vector, shape.x, shape.y);
  run.bind_imagef("field", advected_field.vector, shape.x, shape.y);
  run.bind_imagef("dx", dx.vector, shape.x, shape.y);
  run.bind_imagef("dy", dy.vector, shape.x, shape.y);

  if (p_mask)
  {
    run.bind_imagef("mask", p_mask->vector, shape.x, shape.y);
  }
  else
  {
    mask = Array(shape, 1.f);
    run.bind_imagef("mask", mask.vector, shape.x, shape.y);
  }

  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x, shape.y, advection_length, value_persistence);
  run.execute({shape.x, shape.y});
  run.read_imagef("out");

  return out;
}

Array advection_warp(const Array &z,
                     const Array &advected_field,
                     float        advection_length,
                     float        value_persistence,
                     const Array *p_mask)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, advected_field))
    return Array();
  if (p_mask && !validate_same_shape(z, *p_mask)) return Array();

  Array dx = hmap::gradient_x(z);
  Array dy = hmap::gradient_y(z);

  return advection_warp(z,
                        advected_field,
                        dx,
                        dy,
                        advection_length,
                        value_persistence,
                        p_mask);
}

} // namespace hmap::gpu
