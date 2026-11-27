/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

Array rift(glm::ivec2    shape,
           float         angle,
           float         radius,
           float         axial_slope,
           float         depth,
           bool          scale_with_depth,
           RadialProfile profile,
           float         profile_param,
           float         bottom_extent,
           float         bottom_depth,
           RadialProfile bottom_profile,
           float         bottom_profile_param,
           bool          bottom_force_minimum_depth,
           float         outer_slope,
           const Array  *p_noise_r,
           const Array  *p_noise_offset,
           glm::vec2     center,
           const UvRect &uv_rect,
           Array        *p_rift_mask,
           Array        *p_bottom_mask)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array(shape);
  if (p_noise_offset && !validate_same_shape(shape, *p_noise_offset))
    return Array(shape);

  const float alpha = angle / 180.f * M_PI;
  const float ca = std::cos(alpha);
  const float sa = std::sin(alpha);

  Array z = Array(shape);
  Array rift_mask = Array(shape);
  Array bottom_mask = Array(shape);

  // --- Radial profiles

  auto radial_profile = get_radial_profile_function(profile, profile_param);
  auto bottom_radial_profile = get_radial_profile_function(
      bottom_profile,
      bottom_profile_param);

  // --- Generate

  float bottom_depth_min = 0.f;
  if (bottom_force_minimum_depth) bottom_depth_min = 0.1f * bottom_depth;

  // grid
  bool               endpoint = false;
  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, uv_rect, endpoint);

  // main loop
  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      // --- Compute distance to rift main axis

      float dx = x[i] - center.x;
      float dy = y[j] - center.y;
      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;

      // --- Axial slope

      float da = dx * ca + dy * sa;
      float current_depth = std::max(0.f, depth - axial_slope * da);
      float current_bottom_depth = std::max(bottom_depth_min,
                                            bottom_depth - axial_slope * da);

      float ds = p_noise_offset ? (*p_noise_offset)(i, j) : 0.f;

      // --- Transverse profile

      float current_radius = scale_with_depth
                                 ? radius *
                                       std::max(1e-3f, current_depth / depth)
                                 : radius;

      float d = std::abs(ds + (-dx * sa + dy * ca) / current_radius);
      d = std::max(0.f, d + dr);

      if (d < bottom_extent)
      {
        float dn = d / bottom_extent;
        dn = bottom_radial_profile(dn);
        z(i,
          j) = lerp(-current_depth - current_bottom_depth, -current_depth, dn);

        bottom_mask(i, j) = 1.f - dn;
      }
      else if (d < 1.f)
      {
        float dn = (d - bottom_extent) / (1.f - bottom_extent);
        dn = radial_profile(dn);
        z(i, j) = lerp(-current_depth, 0.f, dn);
      }
      else
      {
        z(i, j) = outer_slope * (d - 1.f) * current_radius;
      }

      rift_mask(i, j) = std::max(0.f, 1.f - d);
    }

  if (p_rift_mask) *p_rift_mask = rift_mask;
  if (p_bottom_mask) *p_bottom_mask = bottom_mask;

  return z;
}

} // namespace hmap
