/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/geo.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

Array basalt_field(glm::ivec2    shape,
                   glm::vec2     kw,
                   std::uint32_t seed,
                   float         warp_kw,
                   float         large_scale_warp_amp,
                   float         large_scale_gain,
                   float         large_scale_amp,
                   float         medium_scale_kw_ratio,
                   float         medium_scale_warp_amp,
                   float         medium_scale_gain,
                   float         medium_scale_amp,
                   float         small_scale_kw_ratio,
                   float         small_scale_amp,
                   float         small_scale_overlay_amp,
                   float         rugosity_kw_ratio,
                   float         rugosity_amp,
                   bool          flatten_activate,
                   float         flatten_kw_ratio,
                   float         flatten_amp,
                   float         bias,
                   const Array  *p_noise_x,
                   const Array  *p_noise_y,
                   glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  // --- large scales

  Array z_large;

  {
    int   octaves = 8;
    float weight = 0.5f;
    float persistence = 0.5f;
    float lacunarity = 2.3f;

    Array dx = gpu::noise_fbm(NoiseType::SIMPLEX2,
                              shape,
                              glm::vec2(warp_kw, warp_kw),
                              seed++,
                              octaves,
                              weight,
                              persistence,
                              lacunarity,
                              nullptr,
                              p_noise_x,
                              p_noise_y,
                              bbox);
    remap(dx, 0.f, large_scale_warp_amp, -1.f, 1.f);

    // base
    glm::vec2         jitter(1.f, 1.f);
    float             k_smoothing = 0.9f;
    VoronoiReturnType return_type = VoronoiReturnType::EDGE_DISTANCE_SQUARED;

    lacunarity = 1.66f;

    z_large = voronoi_fbm(shape,
                          kw,
                          seed++,
                          jitter,
                          bias,
                          k_smoothing,
                          0.f,
                          return_type,
                          octaves,
                          weight,
                          persistence,
                          lacunarity,
                          nullptr,
                          &dx,
                          &dx,
                          bbox);
    remap(z_large, 0.f, 1.f, -0.25, 0.25f);
    z_large = sqrt_safe(z_large);
    gain(z_large, large_scale_gain);

    const float c0 = 0.36f;
    const float c1 = 0.40f;

    z_large = maximum(c1 * z_large, c0);

    // rescale
    remap(z_large, 0.f, large_scale_amp, c0, c1);
  }

  // --- medium scales

  Array z_medium;

  {
    glm::vec2         jitter(1.f, 1.f);
    float             k_smoothing = 0.f;
    VoronoiReturnType return_type = VoronoiReturnType::EDGE_DISTANCE_SQUARED;
    int               octaves = 4;
    float             weight = 0.5f;
    float             persistence = 0.5f;
    float             lacunarity = 2.3f;

    Array dx = gpu::voronoi_fbm(shape,
                                0.5f * glm::vec2(warp_kw, warp_kw),
                                seed++,
                                jitter,
                                /* bias */ 0.f,
                                k_smoothing,
                                0.f,
                                return_type,
                                octaves,
                                weight,
                                persistence,
                                lacunarity,
                                nullptr,
                                p_noise_x,
                                p_noise_y,
                                bbox);
    dx *= medium_scale_warp_amp;

    // base
    lacunarity = 1.7f;
    k_smoothing = 0.9f;

    z_medium = voronoi_fbm(shape,
                           medium_scale_kw_ratio * kw,
                           seed++,
                           jitter,
                           bias,
                           k_smoothing,
                           0.f,
                           return_type,
                           octaves,
                           weight,
                           persistence,
                           lacunarity,
                           nullptr,
                           &dx,
                           &dx,
                           bbox);

    // rescale to [0, 1] (roughly)
    z_medium += 0.25f;
    z_medium *= 2.f;

    z_medium = sqrt_safe(z_medium);
    gain(z_medium, medium_scale_gain);

    remap(z_medium, -1.f + medium_scale_amp, medium_scale_amp, 0.f, 1.f);
  }

  // --- small scales

  Array z_small;

  {
    // some parameters are the same as medium scales

    glm::vec2         jitter(1.f, 1.f);
    float             k_smoothing = 0.f;
    VoronoiReturnType return_type = VoronoiReturnType::EDGE_DISTANCE_SQUARED;
    int               octaves = 4;
    float             weight = 0.5f;
    float             persistence = 0.5f;
    float             lacunarity = 2.f;

    Array dx = gpu::voronoi_fbm(shape,
                                2.f * glm::vec2(warp_kw, warp_kw),
                                seed++,
                                jitter,
                                /* bias */ 0.f,
                                k_smoothing,
                                0.f,
                                return_type,
                                octaves,
                                weight,
                                persistence,
                                lacunarity,
                                nullptr,
                                p_noise_x,
                                p_noise_y,
                                bbox);
    dx *= medium_scale_warp_amp;

    // base
    lacunarity = 1.6f;
    k_smoothing = 0.9f;

    z_small = voronoi_fbm(shape,
                          small_scale_kw_ratio * kw,
                          seed++,
                          jitter,
                          bias,
                          k_smoothing,
                          0.f,
                          return_type,
                          octaves,
                          weight,
                          persistence,
                          lacunarity,
                          nullptr,
                          &dx,
                          &dx,
                          bbox);

    // rescale to [0, 1] (roughly)
    remap(z_small, 0.f, 1.f, -0.25f, 0.25f);

    z_small = sqrt_safe(z_small);
    gain(z_small, medium_scale_gain);
    z_small *= small_scale_amp;
  }

  Array z = hmap::maximum_smooth(z_large, z_medium, 0.04f);

  z *= large_scale_amp / medium_scale_amp;

  z = hmap::maximum_smooth(z, z_small, 0.08f);

  // small scales overlay
  z += small_scale_overlay_amp / small_scale_amp * z_small;

  // --- last steps...

  // multiply add
  {
    int   octaves = 8;
    float weight = 0.5f;
    float persistence = 0.5f;
    float lacunarity = 2.f;

    Array rugosity = gpu::noise_fbm(NoiseType::SIMPLEX2,
                                    shape,
                                    rugosity_kw_ratio * kw,
                                    seed++,
                                    octaves,
                                    weight,
                                    persistence,
                                    lacunarity,
                                    nullptr,
                                    p_noise_x,
                                    p_noise_y,
                                    bbox);
    remap(rugosity, 0.f, 1.f, -1.f, 1.f);

    z += rugosity_amp * rugosity * z;
  }

  // flattening pass
  if (flatten_activate)
  {
    int   octaves = 8;
    float weight = 0.5f;
    float persistence = 0.5f;
    float lacunarity = 2.f;

    Array z_flatten = gpu::noise_fbm(NoiseType::SIMPLEX2,
                                     shape,
                                     flatten_kw_ratio * kw,
                                     seed++,
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity,
                                     nullptr,
                                     p_noise_x,
                                     p_noise_y,
                                     bbox);

    remap(z_flatten, 0.f, 2.f * large_scale_amp + flatten_amp, -1.f, 1.f);
    z = hmap::minimum_smooth(z, z_flatten, 0.3f);
  }

  return z;
}

} // namespace hmap::gpu
