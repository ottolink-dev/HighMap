/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/erosion/coastal_erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

void coastal_erosion_diffusion(Array       &z,
                               Array       &water_depth,
                               float        additional_depth,
                               int          iterations,
                               const Array *p_mask,
                               Array       *p_water_mask)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, water_depth)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  const glm::ivec2 &shape = z.shape;

  Array mask;
  Array z_bckp;

  for (int it = 0; it < iterations; ++it)
  {
    z_bckp = z;
    mask = water_mask(water_depth, z, additional_depth);

    if (p_mask) mask *= (*p_mask);

    // filtering
    hmap::laplace(z, &mask, /* sigma */ 0.125f, 1);

    // adjust water depth so that water height is the same as before
    // filtering
    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        if (water_depth(i, j) > 0.f)
          water_depth(i, j) = z_bckp(i, j) + water_depth(i, j) - z(i, j);
      }
  }

  if (p_water_mask) *p_water_mask = mask;
}

void coastal_erosion_profile(Array       &z,
                             Array       &water_depth,
                             float        shore_ground_extent,
                             float        shore_water_extent,
                             float        slope_shore,
                             float        slope_shore_water,
                             float        scarp_extent_ratio,
                             bool         apply_post_filter,
                             int          post_filter_iterations,
                             bool         solid_shore_mask,
                             float        scarp_mask_transition_ratio,
                             const Array *p_noise,
                             Array       *p_shore_mask,
                             Array       *p_scarp_mask)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, water_depth)) return;
  if (p_noise && !validate_same_shape(z, *p_noise)) return;

  const glm::ivec2 &shape = z.shape;
  Array             z_bckp = z;
  Array             shore_mask(shape); // includes ground & water
  Array             smooth_mask(shape);
  Array             scarp_mask(shape);
  Mat<glm::ivec2>   closest_g(shape); // ground
  Mat<glm::ivec2>   closest_w(shape); // water

  Array r_ground = distance_transform_with_closest(water_depth, closest_g);
  Array r_water = distance_transform_with_closest(is_zero(water_depth),
                                                  closest_w);

  float slope_shore_n = slope_shore / float(shape.x);
  float slope_shore_water_n = slope_shore_water / float(shape.x);
  float t_scarp = 1.f - scarp_extent_ratio;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      if (r_ground(i, j) > 0.f)
      {
        // --- ground

        float dr = p_noise ? (*p_noise)(i, j) : 0.f;

        // transition factor (extent is at least 1 pixel)
        float local_extent = std::max(1.f, shore_ground_extent * (1.f + dr));
        float t = r_ground(i, j) / local_extent;

        if (t <= 1.f)
        {
          shore_mask(i, j) = solid_shore_mask ? 1.f : 1.f - t;
          smooth_mask(i, j) = 1.f - t;

          float zref = z(closest_g(i, j));
          float h = zref + slope_shore_n * r_ground(i, j);

          float new_z = 0.f;

          if (t < t_scarp)
          {
            // shore
            new_z = h;
          }
          else
          {
            // scarp
            float ts = (t - t_scarp) / (1.f - t_scarp); // in [0, 1]
            ts = smoothstep3(ts);

            new_z = lerp(h, z(i, j), ts);

            // sharp mask transition (only a ratio of the width to blend-in)
            scarp_mask(i, j) = threshold_smooth(ts,
                                                0.f,
                                                scarp_mask_transition_ratio);
          }

          z(i, j) = std::min(z(i, j), new_z);
        }
        else if (t <= 1.1f && t_scarp < 1.f)
        {
          // extend the scarp mask slightly outside (10%) to ensure a smooth
          // transition
          float ts = (t - 1.1f) / (1.f - 1.1f); // in [1, 0]
          ts = smoothstep3(ts);
          scarp_mask(i, j) = ts;
        }
      }
      else
      {
        // --- underwater

        // transition factor
        float t = r_water(i, j) / shore_water_extent;

        if (t <= 1.f)
        {
          shore_mask(i, j) = 1.f - t;
          smooth_mask(i, j) = 1.f;

          // ensure slope continuity at water level
          float slope = lerp(slope_shore_n, slope_shore_water_n, t);

          float zref = z(closest_w(i, j));
          float h = zref - slope * r_water(i, j);
          float new_z = lerp(h, z(i, j), smoothstep3(t));

          z(i, j) = std::min(z(i, j), new_z);
        }
      }
    }

  // postprocessing - filter numerical artifacts
  if (apply_post_filter)
  {
    smooth_mask = threshold_smooth(smooth_mask, 1.f - t_scarp, 1.f);
    laplace(z, &smooth_mask, 0.125f, post_filter_iterations);
  }

  // adjust water depth so that water height is the same as before
  // filtering
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      if (water_depth(i, j) > 0.f)
      {
        float water_surface = z_bckp(i, j) + water_depth(i, j);
        water_depth(i, j) = std::max(0.f, water_surface - z(i, j));
      }
      else
      {
        // restore dry cells the inner call may have dirtied
        water_depth(i, j) = 0.f;
      }
    }

  // other optional outputs
  if (p_shore_mask) *p_shore_mask = shore_mask;
  if (p_scarp_mask) *p_scarp_mask = scarp_mask;
}

void coastal_erosion_profile(Array       &z,
                             const Array *p_mask,
                             Array       &water_depth,
                             float        shore_ground_extent,
                             float        shore_water_extent,
                             float        slope_shore,
                             float        slope_shore_water,
                             float        scarp_extent_ratio,
                             bool         apply_post_filter,
                             int          post_filter_iterations,
                             bool         solid_shore_mask,
                             float        scarp_mask_transition_ratio,
                             const Array *p_noise,
                             Array       *p_shore_mask,
                             Array       *p_scarp_mask)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, water_depth)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_noise && !validate_same_shape(z, *p_noise)) return;

  if (!p_mask)
  {
    coastal_erosion_profile(z,
                            water_depth,
                            shore_ground_extent,
                            shore_water_extent,
                            slope_shore,
                            slope_shore_water,
                            scarp_extent_ratio,
                            apply_post_filter,
                            post_filter_iterations,
                            solid_shore_mask,
                            scarp_mask_transition_ratio,
                            p_noise,
                            p_shore_mask,
                            p_scarp_mask);
  }
  else
  {
    Array z_bckp = z;
    Array water_depth_bckp = water_depth;

    Array z_f = z;

    coastal_erosion_profile(z_f,
                            water_depth,
                            shore_ground_extent,
                            shore_water_extent,
                            slope_shore,
                            slope_shore_water,
                            scarp_extent_ratio,
                            apply_post_filter,
                            post_filter_iterations,
                            solid_shore_mask,
                            scarp_mask_transition_ratio,
                            p_noise,
                            p_shore_mask,
                            p_scarp_mask);

    z = lerp(z, z_f, *p_mask);

    // recompute water_depth
    for (int j = 0; j < z.shape.y; ++j)
      for (int i = 0; i < z.shape.x; ++i)
      {
        if (water_depth_bckp(i, j) > 0.f)
        {
          float water_surface = z_bckp(i, j) + water_depth_bckp(i, j);
          water_depth(i, j) = std::max(0.f, water_surface - z(i, j));
        }
        else
        {
          // restore dry cells the inner call may have dirtied
          water_depth(i, j) = 0.f;
        }
      }

    // output masks
    if (p_shore_mask) *p_shore_mask = (*p_shore_mask) * (*p_mask);
    if (p_scarp_mask) *p_scarp_mask = (*p_scarp_mask) * (*p_mask);
  }
}

} // namespace hmap
