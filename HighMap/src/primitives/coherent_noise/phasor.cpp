/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>
#include <functional>

#include "highmap/array.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

Array phasor(PhasorProfile   phasor_profile,
             glm::ivec2      shape,
             float           kp_global,
             std::uint32_t   seed,
             float           angle_shift,
             float           normalization,
             const glm::vec2 jitter,
             float           delta,
             float           phase_smoothing,
             const Array    *p_angle,
             const Array    *p_noise_x,
             const Array    *p_noise_y,
             glm::vec4       bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_angle && !validate_same_shape(shape, *p_angle)) return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  // wavenumbers
  float           kp = 1.f;
  const glm::vec2 kw = {kp_global, kp_global};

  // angle field
  Array angle(shape, angle_shift);
  if (p_angle) angle += *p_angle;

  // compute phase
  Array modulus(shape);
  Array phase = phase_field_angle(angle,
                                  kw,
                                  seed,
                                  kp,
                                  normalization,
                                  jitter,
                                  /* p_ctrl_param */ nullptr,
                                  p_noise_x,
                                  p_noise_y,
                                  &modulus,
                                  /* p_angle_jump_mask */ nullptr,
                                  bbox);

  // apply phase function
  float profile_avg = 0.f;
  auto  fct = get_phasor_profile_function(phasor_profile, delta, &profile_avg);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      float value = fct(phase(i, j));

      // modulus-based smoothing to avoid kinky transitions between
      // phases
      float t = 2.f / M_PI * std::atan(phase_smoothing * modulus(i, j));

      phase(i, j) = lerp(profile_avg, value, t);
    }

  return phase;
}

Array phasor_fbm(PhasorProfile   phasor_profile,
                 glm::ivec2      shape,
                 float           kp_global,
                 std::uint32_t   seed,
                 float           angle_shift,
                 int             octaves,
                 float           weight,
                 float           persistence,
                 float           lacunarity,
                 float           normalization,
                 const glm::vec2 jitter,
                 float           delta,
                 float           phase_smoothing,
                 const Array    *p_angle,
                 const Array    *p_noise_x,
                 const Array    *p_noise_y,
                 glm::vec4       bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_angle && !validate_same_shape(shape, *p_angle)) return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array sum(shape, 0.f);
  Array amp(shape, 1.f);

  for (int k = 0; k < octaves; k++)
  {
    Array value = phasor(phasor_profile,
                         shape,
                         kp_global,
                         seed,
                         angle_shift,
                         normalization,
                         jitter,
                         delta,
                         phase_smoothing,
                         p_angle,
                         p_noise_x,
                         p_noise_y,
                         bbox);

    sum += value * amp;
    amp *= (1.f - weight) + weight * minimum(value + 1.f, 2.f) * 0.5f;
    kp_global *= lacunarity;
    amp *= persistence;
    seed++;
  }

  return sum;
}

} // namespace hmap::gpu
