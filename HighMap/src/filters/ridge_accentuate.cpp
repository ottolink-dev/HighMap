/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

void ridge_accentuate(Array &array, float strength, int ir, bool reverse)
{
  if (!validate_non_empty(array)) return;

  Array array_f = array;
  Array array_new = Array(array.shape);

  // --- Smooth terrain for normal computation

  if (ir > 0) smooth_cpulse(array_f, ir);

  if (reverse) strength = -strength;

  int nx = array.shape.x;
  int ny = array.shape.y;

  float scale_x = (float)(nx - 1);
  float scale_y = (float)(ny - 1);

  // --- Displaced sampling along terrain normal direction

  for (int j = 0; j < ny; j++)
    for (int i = 0; i < nx; i++)
    {
      // central difference scaled to domain units
      float dx = 0.5f *
                 (array_f(std::min(i + 1, nx - 1), j) -
                  array_f(std::max(i - 1, 0), j)) *
                 scale_x;
      float dy = 0.5f *
                 (array_f(i, std::min(j + 1, ny - 1)) -
                  array_f(i, std::max(j - 1, 0))) *
                 scale_y;

      glm::vec3 n(-dx, -dy, 1.f);
      n /= std::hypot(n.x, n.y, n.z);

      // displacement along horizontal normal direction
      // sampling position in pixel coordinates
      float x = (float)i + strength * scale_x * n.x;
      float y = (float)j + strength * scale_y * n.y;

      // clamp sampling coordinates to array bounds
      x = std::clamp(x, 0.f, scale_x);
      y = std::clamp(y, 0.f, scale_y);

      // bilinear interpolation parameters
      int ip = std::min((int)x, nx - 2);
      int jp = std::min((int)y, ny - 2);
      if (ip < 0) ip = 0;
      if (jp < 0) jp = 0;

      float u = x - (float)ip;
      float v = y - (float)jp;

      array_new(i, j) = array.get_value_bilinear_at(ip, jp, u, v);
    }

  array = array_new;
}

void ridge_accentuate(Array       &array,
                      const Array *p_mask,
                      float        strength,
                      int          ir,
                      bool         reverse)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { ridge_accentuate(a, strength, ir, reverse); });
}

} // namespace hmap
