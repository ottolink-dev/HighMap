/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap::gpu
{

Array coastal_fetch(const Array &z,
                    int          ndirections,
                    const Array *p_compute_mask)
{
  if (!validate_non_empty(z)) return Array();
  if (p_compute_mask && !validate_same_shape(z, *p_compute_mask))
    return Array();

  const glm::ivec2 &shape = z.shape;
  Array             out(shape);
  Array             compute_mask(shape, 1.f);

  if (p_compute_mask) compute_mask = *p_compute_mask;

  auto run = clwrapper::Run("coastal_fetch");

  run.bind_imagef("z", z.vector, shape.x, shape.y);
  run.bind_imagef("compute_mask", compute_mask.vector, shape.x, shape.y);
  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x, shape.y, ndirections);

  run.execute({shape.x, shape.y});

  run.read_imagef("out");

  return out;
}

Array coastal_fetch_directional(const Array &z,
                                float        angle,
                                float        directional_exp,
                                int          ndirections,
                                const Array *p_compute_mask)
{
  if (!validate_non_empty(z)) return Array();
  if (p_compute_mask && !validate_same_shape(z, *p_compute_mask))
    return Array();

  const glm::ivec2 &shape = z.shape;
  Array             out(shape);
  Array             compute_mask(shape, 1.f);

  float                alpha = angle / 180.f * M_PI;
  std::array<float, 2> wind_dir = {std::cos(alpha), std::sin(alpha)};

  if (p_compute_mask) compute_mask = *p_compute_mask;

  auto run = clwrapper::Run("coastal_fetch_directional");

  run.bind_imagef("z", z.vector, shape.x, shape.y);
  run.bind_imagef("compute_mask", compute_mask.vector, shape.x, shape.y);
  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x, shape.y, ndirections, wind_dir, directional_exp);

  run.execute({shape.x, shape.y});

  run.read_imagef("out");

  return out;
}

} // namespace hmap::gpu
