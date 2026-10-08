/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/primitives/geo.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

Array mountain_inselberg(glm::ivec2    shape,
                         std::uint32_t seed,
                         float         scale,
                         int           octaves,
                         float         rugosity,
                         float         angle,
                         float         gamma,
                         bool          round_shape,
                         bool          add_deposition,
                         float         bulk_amp,
                         float         base_noise_amp,
                         float         bias,
                         float         k_smoothing,
                         glm::vec2     center,
                         const Array  *p_noise_x,
                         const Array  *p_noise_y,
                         glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  // apply global scaling to reference values
  const float     half_width = 0.2f * scale;
  const glm::vec2 kw = glm::vec2(2.6f / scale, 2.6f / scale);

  const float persistence = 0.5f;
  const float lacunarity = 2.f;
  const float alpha = angle / 180.f * M_PI;

  // prepare base noise used for displacements
  Array noise = scale * base_noise_amp *
                gpu::noise_fbm(NoiseType::SIMPLEX2,
                               shape,
                               kw,
                               seed,
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

  // enveloppe pulse
  Array *p_gx = round_shape ? nullptr : &dx;
  Array *p_gy = round_shape ? nullptr : &dy;

  Array pulse = gaussian_pulse(shape,
                               half_width,
                               /* p_ctrl_param */ nullptr,
                               p_gx,
                               p_gy,
                               center,
                               bbox);

  // base primitives
  glm::vec2         jitter(1.f, 1.f);
  VoronoiReturnType return_type = VoronoiReturnType::CONSTANT_F2MF1_SQUARED;

  Array voronoi = 0.72f + voronoi_fbm(shape,
                                      kw,
                                      seed,
                                      jitter,
                                      bias,
                                      k_smoothing,
                                      0.f,
                                      return_type,
                                      octaves,
                                      /* weight */ 0.7f,
                                      persistence,
                                      lacunarity,
                                      nullptr,
                                      &dx,
                                      &dy,
                                      bbox);

  clamp_min(voronoi, 0.f);

  voronoi *= pulse;

  if (bulk_amp > 0.f)
  {
    voronoi += bulk_amp * pulse;
    voronoi *= 1.f / (1.f + bulk_amp);
  }

  gamma_correction(voronoi, gamma);

  if (add_deposition)
  {
    int   ir = (int)(0.05f * scale * shape.x);
    float k = 0.05f;
    gpu::smooth_fill(voronoi, ir, k);
  }

  return voronoi;
}

} // namespace hmap::gpu
