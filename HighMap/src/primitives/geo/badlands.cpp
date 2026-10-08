/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/geo.hpp"

namespace hmap::gpu
{

Array badlands(glm::ivec2    shape,
               glm::vec2     kw,
               std::uint32_t seed,
               int           octaves,
               float         rugosity,
               float         angle,
               float         bias,
               float         k_smoothing,
               float         base_noise_amp,
               const Array  *p_noise_x,
               const Array  *p_noise_y,
               glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  const float persistence = 0.5f;
  const float lacunarity = 2.3f;
  const float alpha = angle / 180.f * M_PI;

  // prepare base noise used for displacements
  Array noise = base_noise_amp * gpu::noise_fbm(NoiseType::SIMPLEX2,
                                                shape,
                                                kw,
                                                seed++,
                                                octaves,
                                                rugosity,
                                                persistence,
                                                lacunarity,
                                                nullptr,
                                                p_noise_x,
                                                p_noise_y,
                                                bbox);

  Array dx = noise * std::cos(alpha);
  Array dy = noise * std::sin(alpha);

  // base primitive
  glm::vec2         jitter(1.f, 1.f);
  VoronoiReturnType return_type = VoronoiReturnType::CONSTANT_F2MF1_SQUARED;

  Array voronoi = voronoi_fbm(shape,
                              kw,
                              seed++,
                              jitter,
                              bias,
                              k_smoothing,
                              0.f,
                              return_type,
                              octaves,
                              /* weight */ 0.5f,
                              persistence,
                              lacunarity,
                              nullptr,
                              &dx,
                              &dy,
                              bbox);

  return voronoi;
}

} // namespace hmap::gpu
