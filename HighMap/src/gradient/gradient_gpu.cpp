/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

Array gradient_angle_circular_smoothing(const Array &array,
                                        int          ir,
                                        bool         downward)
{
  if (!validate_non_empty(array)) return Array();

  // gradients
  Array dx = gradient_x(array);
  Array dy = gradient_y(array);
  Array dn = hypot(dx, dy);
  Array dn_safe = maximum(dn, 1e-9f);

  if (downward)
  {
    dx = -1.f * dx;
    dy = -1.f * dy;
  }

  // angle to unit vector
  Array u = dx / dn_safe;
  Array v = dy / dn_safe;

  // smoothing
  gpu::smooth_cpulse(u, ir);
  gpu::smooth_cpulse(v, ir);
  gpu::smooth_cpulse(dn, ir);

  // renormalize and compute the angle
  dn_safe = maximum(dn, 1e-9f);
  u /= dn_safe;
  v /= dn_safe;

  return atan2(v, u);
}

Array laplacian_fract(const Array &array, float s, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array out(array.shape);

  auto run = clwrapper::Run("laplacian_fract");

  run.bind_imagef("array",
                  const_cast<std::vector<float> &>(array.vector),
                  array.shape.x,
                  array.shape.y);
  run.bind_imagef("out", out.vector, array.shape.x, array.shape.y, true);

  run.bind_arguments(array.shape.x, array.shape.y, ir, s);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return out;
}

} // namespace hmap::gpu
