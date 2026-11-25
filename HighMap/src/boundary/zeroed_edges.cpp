/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <functional>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/distance_functions.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

void zeroed_edges(Array               &array,
                  RadialProfile        radial_profile,
                  float                profile_param,
                  float                amount,
                  DistanceFunction     distance,
                  DistanceFunctionAxis dfa,
                  glm::vec2            center,
                  float                radius,
                  const Array         *p_noise_r,
                  const UvRect        &uv_rect)
{
  if (!validate_non_empty(array)) return;
  if (p_noise_r && !validate_same_shape(array, *p_noise_r)) return;
  if (!validate_not_zero(radius)) return;

  const glm::ivec2 &shape = array.shape;

  auto radial_fct = get_radial_profile_function(radial_profile, profile_param);
  auto dist_fct = get_distance_function(distance, dfa);

  // grid
  bool               endpoint = false;
  std::vector<float> x, y;
  grid_xy_vector(x, y, shape, uv_rect, endpoint);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      // --- Compute radius

      float dx = x[i] - center.x;
      float dy = y[j] - center.y;
      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
      float r = dr + dist_fct(dx, dy) / radius;
      r = std::clamp(r, 0.f, 1.f);

      // --- Apply envelope

      array(i, j) *= lerp(1.f, 1.f - radial_fct(r), amount);
    }
}

} // namespace hmap
