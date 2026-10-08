/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/morphology.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/geo.hpp"
#include "highmap/range.hpp"

namespace hmap
{

Array island_land_mask(glm::ivec2       shape,
                       float            radius,
                       std::uint32_t    seed,
                       float            displacement,
                       NoiseType        noise_type,
                       float            kw,
                       int              octaves,
                       float            weight,
                       float            persistence,
                       float            lacunarity,
                       const glm::vec2 &center,
                       const glm::vec4 &bbox)
{
  if (!validate_shape(shape)) return Array();

  Array mask(shape);

  std::unique_ptr<NoiseFunction> p = create_noise_function_from_type(
      noise_type,
      glm::vec2(0.5f * kw, 0.5f * kw),
      seed);

  FbmFunction f = FbmFunction(std::move(p),
                              octaves,
                              weight,
                              persistence,
                              lacunarity);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      float x = float(i) / float(shape.x - 1);
      float y = float(j) / float(shape.y - 1);

      x = (bbox.y - bbox.x) * x + bbox.x - center.x;
      y = (bbox.w - bbox.z) * y + bbox.z - center.y;

      float r = std::hypot(x, y);
      float theta = std::atan2(y, x);

      float dr = displacement *
                 f.get_delegate()(std::cos(theta), std::sin(theta), 0.f);

      mask(i, j) = r < radius + dr ? 1.f : 0.f;
    }

  return mask;
}

} // namespace hmap

namespace hmap::gpu
{

inline float helper_radial_profile(float r,
                                   float slope_start,
                                   float slope_end,
                                   float apex_elevation,
                                   float k_smooth,
                                   float radial_gain)
{
  float hs = slope_start * r;
  float he = hmap::maximum_smooth(0.f,
                                  slope_end * (r - 1.f) + apex_elevation,
                                  k_smooth);
  return lerp(hs, he, gain(r, radial_gain));
}

inline float helper_apply_leeward(float h,
                                  float r,
                                  float hs,
                                  float slope_max,
                                  float k_smooth,
                                  float lee_amp,
                                  float alpha,
                                  float lee_alpha)
{
  float hl = h + lee_amp * smoothstep3(1.f - r) * std::cos(alpha - lee_alpha);
  hl = hmap::maximum_smooth(hl, hs, k_smooth);
  hl = hmap::minimum_smooth(hl, slope_max * r, k_smooth);
  return std::max(0.f, hl);
}

inline float helper_apply_uplift(float h,
                                 float r,
                                 float slope_max,
                                 float uplift_amp,
                                 float k_smooth)
{
  float hl = h + uplift_amp;
  hl = hmap::minimum_smooth(hl, slope_max * r, k_smooth);
  return std::max(0.f, hl);
}

Array island(const Array &land_mask,
             const Array *p_noise_r,
             float        apex_elevation,
             bool         filter_distance,
             int          filter_ir,
             float        slope_min,
             float        slope_max,
             float        slope_start,
             float        slope_end,
             float        slope_noise_intensity,
             float        k_smooth,
             float        radial_noise_intensity,
             float        helper_radial_profile_gain,
             float        water_decay,
             float        water_depth,
             float        lee_angle,
             float        lee_amp,
             float        uplift_amp,
             Array       *p_water_depth,
             Array       *p_inland_mask)
{
  if (!validate_non_empty(land_mask)) return Array();
  if (p_noise_r && !validate_same_shape(land_mask, *p_noise_r))
    return Array(land_mask.shape);

  glm::ivec2 shape = land_mask.shape;
  float      lee_alpha = lee_angle * (float(M_PI) / 180.f);

  // --- distance from land
  Array r_ground = distance_transform(is_zero(land_mask));
  if (filter_distance) gpu::smooth_cpulse(r_ground, filter_ir);

  float rmax;
  int   ic, jc;
  r_ground.argmax(rmax, ic, jc);

  // --- PASS 1: normalize radius + apply radial noise
  const bool use_noise = p_noise_r != nullptr;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float r = r_ground(i, j) / rmax;

      if (use_noise)
      {
        float dr = radial_noise_intensity * (*p_noise_r)(i, j);
        r = std::max(0.f, r - std::max(0.f, dr));
      }

      r_ground(i, j) = std::min(r, 1.f);
    }

  // --- water distance
  Array r_water = distance_transform(r_ground);

  // scale water dist
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      r_water(i, j) /= (water_decay * float(shape.x));

  // --- prepare outputs
  Array z_ground(shape);
  Array z_water(shape);

  // --- PASS 2: compute ground + underwater topography
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float r = r_ground(i, j);
      float r_w = r_water(i, j);

      z_water(i, j) = std::min(0.f, -water_depth * (1.f - std::exp(-r_w)));

      if (r == 0.f)
      {
        z_ground(i, j) = 0.f;
        continue;
      }

      // slope noise
      float slope0 = slope_start;
      if (use_noise)
      {
        float dslope = slope_noise_intensity *
                       std::max(0.f, (*p_noise_r)(i, j)) * (1.f - r);

        slope0 = hmap::maximum_smooth(slope_min,
                                      slope0 * (1.f - dslope),
                                      k_smooth);
      }

      float h = helper_radial_profile(r,
                                      slope0,
                                      slope_end,
                                      apex_elevation,
                                      k_smooth,
                                      helper_radial_profile_gain);

      // leeward
      if (lee_amp > 0.f)
      {
        float alpha = std::atan2(j - jc, i - ic);
        h = helper_apply_leeward(h,
                                 r,
                                 slope0 * r,
                                 slope_max,
                                 k_smooth,
                                 lee_amp,
                                 alpha,
                                 lee_alpha);
      }

      // uplift
      if (uplift_amp > 0.f)
        h = helper_apply_uplift(h, r, slope_max, uplift_amp, k_smooth);

      z_ground(i, j) = h;
    }

  if (p_water_depth) *p_water_depth = -z_water;
  if (p_inland_mask) *p_inland_mask = r_ground;

  return z_ground + z_water;
}

Array island(const Array  &land_mask,
             std::uint32_t seed,
             float         noise_amp,
             glm::vec2     noise_kw,
             int           noise_octaves,
             float         noise_rugosity,
             float         noise_angle,
             float         noise_k_smoothing,
             float         apex_elevation,
             bool          filter_distance,
             int           filter_ir,
             float         slope_min,
             float         slope_max,
             float         slope_start,
             float         slope_end,
             float         slope_noise_intensity,
             float         k_smooth,
             float         radial_noise_intensity,
             float         helper_radial_profile_gain,
             float         water_decay,
             float         water_depth,
             float         lee_angle,
             float         lee_amp,
             float         uplift_amp,
             Array        *p_water_depth,
             Array        *p_inland_mask)
{
  if (!validate_non_empty(land_mask)) return Array();

  glm::ivec2 shape = land_mask.shape;

  // --- generate ready-made noise

  hmap::Array base = noise_amp * hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                                 shape,
                                                 noise_kw,
                                                 seed++,
                                                 noise_octaves,
                                                 0.f);

  const float alpha = noise_angle / 180.f * M_PI;
  Array       dx = base * std::cos(alpha);
  Array       dy = base * std::sin(alpha);

  auto noise = hmap::gpu::voronoi_fbm(shape,
                                      noise_kw,
                                      seed++,
                                      {1.f, 1.f},
                                      /* bias */ 0.f,
                                      noise_k_smoothing,
                                      0.f,
                                      hmap::VoronoiReturnType::F2MF1_SQUARED,
                                      noise_octaves,
                                      /* weight */ noise_rugosity,
                                      0.5f,
                                      2.f,
                                      /* p_ctrl_param */ nullptr,
                                      &dx,
                                      &dy);

  remap(noise);

  // --- apply

  return island(land_mask,
                &noise,
                apex_elevation,
                filter_distance,
                filter_ir,
                slope_min,
                slope_max,
                slope_start,
                slope_end,
                slope_noise_intensity,
                k_smooth,
                radial_noise_intensity,
                helper_radial_profile_gain,
                water_decay,
                water_depth,
                lee_angle,
                lee_amp,
                uplift_amp,
                p_water_depth,
                p_inland_mask);
}

} // namespace hmap::gpu
