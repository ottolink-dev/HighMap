/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include "highmap/array.hpp"
#include "highmap/erosion/deprecated.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

void hydraulic_diffusion(Array &z,
                         float  c_diffusion,
                         float  talus,
                         int    iterations)
{
  if (!validate_non_empty(z)) return;

  Array dx = Array(z.shape);
  Array dy = Array(z.shape);
  Array qx = Array(z.shape);
  Array qy = Array(z.shape);
  Array c = Array(z.shape);

  for (int it = 0; it < iterations; it++)
  {
    gradient_x(z, dx);
    gradient_y(z, dy);

    c = 1.f / (1.f - dx * dx / (talus * talus));
    qx = c_diffusion * c * dx;

    c = 1.f / (1.f - dy * dy / (talus * talus));
    qy = c_diffusion * c * dy;

    gradient_x(qx, dx);
    gradient_y(qy, dy);

    z += dx + dy;
  }
}

} // namespace hmap
