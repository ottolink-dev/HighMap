/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/functions.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

void warp(Array &array, const Array *p_dx, const Array *p_dy)
{
  if (!validate_non_empty(array)) return;
  if (p_dx && !validate_same_shape(array, *p_dx)) return;
  if (p_dy && !validate_same_shape(array, *p_dy)) return;

  hmap::ArrayFunction f = hmap::ArrayFunction(array, glm::vec2(1.f, 1.f), true);

  fill_array_using_xy_function(array,
                               {0.f, 1.f, 0.f, 1.f},
                               nullptr,
                               p_dx,
                               p_dy,
                               nullptr,
                               f);
}

void warp_directional(Array &array,
                      float  angle,
                      float  amount,
                      int    ir,
                      bool   reverse)
{
  if (!validate_non_empty(array)) return;

  // --- define displacement: same as warp_downslope but add a
  // --- reference angle to scale the amount of warping
  float angle_rad = angle / 180.f * M_PI;

  Array array_new = Array(array.shape);
  Array alpha = Array(array.shape);

  {
    Array array_f = array;
    if (ir > 0) smooth_cpulse(array_f, ir);
    alpha = gradient_angle(array_f);
  }

  if (reverse) amount = -amount;

  Array ca = amount * cos(alpha - angle_rad);
  Array sa = amount * sin(alpha - angle_rad);

  // --- warp
  warp(array, &ca, &sa);
}

void warp_directional(Array       &array,
                      float        angle,
                      const Array *p_mask,
                      float        amount,
                      int          ir,
                      bool         reverse)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  if (!p_mask)
    warp_directional(array, angle, amount, ir, reverse);
  else
  {
    Array array_f = array;
    warp_directional(array_f, angle, amount, ir, reverse);
    array = lerp(array, array_f, *(p_mask));
  }
}

void warp_downslope(Array &array, float amount, int ir, bool reverse)
{
  if (!validate_non_empty(array)) return;

  Array array_new = Array(array.shape);
  Array alpha = Array(array.shape);

  {
    Array array_f = array;
    if (ir > 0) smooth_cpulse(array_f, ir);

    alpha = gradient_angle(array_f);
  }

  if (reverse) amount = -amount;

  Array ca = amount * cos(alpha);
  Array sa = amount * sin(alpha);

  warp(array, &ca, &sa);
}

void warp_downslope(Array       &array,
                    const Array *p_mask,
                    float        amount,
                    int          ir,
                    bool         reverse)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  if (!p_mask)
    warp_downslope(array, amount, ir, reverse);
  else
  {
    Array array_f = array;
    warp_downslope(array_f, amount, ir, reverse);
    array = lerp(array, array_f, *(p_mask));
  }
}

} // namespace hmap
