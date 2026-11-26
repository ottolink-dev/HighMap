/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <functional>

#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

Array flatbed_carve(glm::ivec2    shape,
                    const Path   &path,
                    float         bottom_extent,
                    float         vmin,
                    float         depth,
                    float         falloff_distance,
                    float         outer_slope,
                    bool          preserve_bedshape,
                    RadialProfile radial_profile,
                    float         radial_profile_parameter,
                    Array        *p_falloff_mask,
                    const Array  *p_noise_r,
                    const UvRect &uv_rect)
{
  if (!validate_shape(shape)) return Array();

  // out
  Array z(shape);
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return z;
  if (!validate_non_empty(path, "Path")) return z;

  Array mask(shape);

  // project path to the array
  Path path_cpy = path;
  path_cpy.set_values(1.f);
  path_cpy.to_array(z, BBox2D::from_vec4(uv_rect.to_vec4()));

  // radial profile
  auto profile_fct = get_radial_profile_function(radial_profile,
                                                 radial_profile_parameter);

  // --- distance and distance function

  Array dist = distance_transform(z);

  std::function<float(int i, int j)> radius_fct;

  if (preserve_bedshape)
  {
    radius_fct = [=](int i, int j)
    {
      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
      float r = dist(i, j) + bottom_extent * dr * (dist(i, j) / bottom_extent);
      return std::max(0.f, r);
    };
  }
  else
  {
    radius_fct = [=](int i, int j)
    {
      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
      float r = dist(i, j) + bottom_extent * dr;
      return std::max(0.f, r);
    };
  }

  // --- actual shape

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float r = radius_fct(i, j);

      if (r < bottom_extent)
      {
        float t = r / bottom_extent; // in [0, 1]
        z(i, j) = vmin + depth * profile_fct(t);
      }
      else
      {
        z(i, j) = vmin + depth +
                  outer_slope * (r - bottom_extent) / float(shape.x - 1);
      }
    }

  // --- generate mask

  if (p_falloff_mask)
  {
    *p_falloff_mask = Array(shape);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        float r = radius_fct(i, j);

        if (r < bottom_extent)
        {
          (*p_falloff_mask)(i, j) = 1.f;
        }
        else if (r < bottom_extent + falloff_distance)
        {
          float t = 1.f - (r - bottom_extent) / falloff_distance;
          (*p_falloff_mask)(i, j) = smoothstep3(t);
        }
      }
  }

  return z;
}

void flatbed_carve(Array        &z,
                   const Path   &path,
                   float         bottom_extent,
                   float         vmin,
                   float         depth,
                   float         falloff_distance,
                   float         outer_slope,
                   bool          preserve_bedshape,
                   RadialProfile radial_profile,
                   float         radial_profile_parameter,
                   Array        *p_falloff_mask,
                   const Array  *p_noise_r,
                   const UvRect &uv_rect)
{
  if (!validate_non_empty(z)) return;
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return;
  if (!validate_non_empty(path, "Path")) return;

  hmap::Array mask;

  Array flatbed = flatbed_carve(z.shape,
                                path,
                                bottom_extent,
                                vmin,
                                depth,
                                falloff_distance,
                                outer_slope,
                                preserve_bedshape,
                                radial_profile,
                                radial_profile_parameter,
                                &mask,
                                p_noise_r,
                                uv_rect);

  z = lerp(z, flatbed, mask);

  if (p_falloff_mask) *p_falloff_mask = mask;
}

} // namespace hmap
