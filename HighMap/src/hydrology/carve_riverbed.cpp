/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/morphology.hpp"
#include "highmap/transform.hpp"

namespace hmap
{

Array carve_riverbed(const Array  &z,
                     const Array  &z_river,
                     float         talus_riverbank,
                     bool          smooth_river_bottom,
                     float         merging_distance,
                     std::uint32_t seed,
                     float         riverbank_noise_ratio,
                     const Array  *p_noise_x,
                     const Array  *p_noise_y)
{
  if (!validate_non_empty(z)) return Array();
  if (!validate_same_shape(z, z_river)) return Array();
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return Array();

  const glm::ivec2 shape = z.shape;
  Array            mask(shape);
  Array            zr = z_river;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if (z_river(i, j) != z(i, j)) mask(i, j) = 1.f;

  const int ir = 2;
  expand_talus(zr, mask, talus_riverbank, seed, ir, riverbank_noise_ratio);

  if (smooth_river_bottom) laplace(zr);

  // Transition mask
  mask = distance_transform(mask);
  mask = exp(-mask / merging_distance);
  laplace(mask);

  // Avoid boundary numerical artefacts
  extrapolate_borders(mask, 2 * ir);

  warp(mask, p_noise_x, p_noise_y);

  return lerp(z, zr, mask);
}

} // namespace hmap
