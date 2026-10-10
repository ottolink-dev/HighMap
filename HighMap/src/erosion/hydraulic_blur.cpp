/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>

#include "highmap/array.hpp"
#include "highmap/erosion/deprecated.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void hydraulic_blur(Array &z, float vmax, float radius, float k_smoothing)
{
  if (!validate_non_empty(z)) return;

  int ir = std::max(1, (int)(radius * z.shape.x));

  Array z_clamped = z;
  if (k_smoothing)
    clamp_min_smooth(z_clamped, vmax, k_smoothing);
  else
    clamp_min_smooth(z_clamped, vmax);

  Array z_filtered = z;
  smooth_cpulse(z_filtered, ir);

  z = (z_clamped + z_filtered - vmax);
}

} // namespace hmap
