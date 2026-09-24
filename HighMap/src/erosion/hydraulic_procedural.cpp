/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>
#include <functional>

#include "highmap/array.hpp"
#include "highmap/blending.hpp"
#include "highmap/erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/functions.hpp"
#include "highmap/gradient.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/morphology.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

void hydraulic_procedural(Array         &z,
                          float          kp_global,
                          float          c_erosion,
                          std::uint32_t  seed,
                          ErosionProfile erosion_profile,
                          float          erosion_profile_parameter,
                          float          angle_shift,
                          float          phase_smoothing,
                          float          talus_ref,
                          float          gradient_scaling_ratio,
                          float          gradient_power,
                          bool           exclude_ridges,
                          bool           apply_deposition,
                          float          deposition_strength,
                          bool           enable_default_noise,
                          float          noise_amp,
                          const Array   *p_kp_multiplier,
                          const Array   *p_angle_shift,
                          const Array   *p_noise_x,
                          const Array   *p_noise_y,
                          Array         *p_ridge_mask,
                          glm::vec4      bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_kp_multiplier && !validate_same_shape(z, *p_kp_multiplier)) return;
  if (p_angle_shift && !validate_same_shape(z, *p_angle_shift)) return;
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return;
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return;

  // ---Resolve derived parameters

  const glm::ivec2 shape = z.shape;
  const int        kp_ir = int(shape.x / kp_global);
  const int        angle_filter_ir = kp_ir;
  const bool       rotate90 = false;
  const float      normalization = 1.f;
  const glm::vec2  jitter = {1.f, 1.f};
  const int        gradient_prefilter_ir = kp_ir;

  // --- Default noise

  Array dx, dy;

  if (enable_default_noise)
  {
    const glm::vec2 kw_noise = {kp_global, kp_global};

    dx = noise_amp * gpu::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                    z.shape,
                                    kw_noise,
                                    ++seed,
                                    /* octaves */ 8,
                                    /* weight*/ 0.7f,
                                    /* persistence */ 0.5f,
                                    /* lacunarity */ 2.f,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    bbox);

    dy = noise_amp * gpu::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                    z.shape,
                                    kw_noise,
                                    ++seed,
                                    /* octaves */ 8,
                                    /* weight*/ 0.7f,
                                    /* persistence */ 0.5f,
                                    /* lacunarity */ 2.f,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    bbox);

    p_noise_x = &dx;
    p_noise_y = &dy;
  }

  // --- Compute phase field

  Array modulus(shape);

  Array phase = gpu::phase_field(z,
                                 seed,
                                 kp_global,
                                 rotate90,
                                 normalization,
                                 jitter,
                                 angle_filter_ir,
                                 // nullptr,
                                 p_kp_multiplier,
                                 p_noise_x,
                                 p_noise_y,
                                 &modulus,
                                 /* p_angle_jump_mask */ nullptr,
                                 bbox);

  // add optional local or global angle shifts
  if (angle_shift != 0.) phase += angle_shift / 180.f * M_PI;

  if (p_angle_shift) phase += *p_angle_shift;

  // --- Ridge profile

  float profile_avg = 0.f;
  auto  pfct = get_erosion_profile_function(erosion_profile,
                                           erosion_profile_parameter,
                                           profile_avg);

  Array ridge(shape);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      float value = pfct(phase(i, j));

      // modulus-based smoothing to avoid kinky transitions between
      // phases
      float t = 2.f / M_PI * std::atan(phase_smoothing * modulus(i, j));

      ridge(i, j) = lerp(profile_avg, value, t);
    }

  // --- Erosion amplitude

  Array zb = z;
  // zb = hmap::bulkify(zb, hmap::PrimitiveType::PRIM_SMOOTH_COSINE, 1.f);

  Array amp = flow_accumulation_dinf(zb, talus_ref);
  amp = log10(amp);

  // scale erosion with local gradient
  Array gn = hmap::gradient_norm(zb);

  {
    gpu::smooth_cpulse(gn, gradient_prefilter_ir);
    remap(gn);
    gn = pow(gn, gradient_power);
    gn = smoothstep5_lower(gn);
    amp *= (1.f - gradient_scaling_ratio) + gradient_scaling_ratio * gn;
  }

  gpu::smooth_cpulse(amp, kp_ir);
  remap(amp);

  ridge *= amp;

  // --- remove line ridges and sinks to avoid artifacts

  if (exclude_ridges)
  {
    Array mr = 1.f - gpu::morphological_top_hat(z, kp_ir); // ridges
    remap(mr);
    ridge *= mr;
  }

  // --- erode

  z -= c_erosion * ridge;

  // --- mimic deposition

  if (apply_deposition)
  {
    int   deposition_ir = kp_ir;
    Array zd = z;
    gpu::smooth_fill_holes(zd, deposition_ir);
    zd = gpu::blend_gradients(zd, z, deposition_ir);
    z = lerp(z, zd, deposition_strength);
  }

  // --- optional output

  if (p_ridge_mask) *p_ridge_mask = cos(phase);
}

void hydraulic_procedural_fbm(Array         &z,
                              float          kp_global,
                              float          c_erosion,
                              std::uint32_t  seed,
                              ErosionProfile erosion_profile,
                              int            octaves,
                              float          persistence,
                              float          lacunarity,
                              float          erosion_profile_parameter,
                              float          angle_shift,
                              float          phase_smoothing,
                              float          talus_ref,
                              float          gradient_scaling_ratio,
                              float          gradient_power,
                              bool           exclude_ridges,
                              bool           apply_deposition,
                              float          deposition_strength,
                              bool           enable_default_noise,
                              float          noise_amp,
                              const Array   *p_kp_multiplier,
                              const Array   *p_angle_shift,
                              const Array   *p_noise_x,
                              const Array   *p_noise_y,
                              Array         *p_ridge_mask,
                              glm::vec4      bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_kp_multiplier && !validate_same_shape(z, *p_kp_multiplier)) return;
  if (p_angle_shift && !validate_same_shape(z, *p_angle_shift)) return;
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return;
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return;

  float a = 1.f;
  float m = 1.f;

  Array ridge_mask = Array(z.shape);
  float a_sum = 0.f;

  for (int k = 0; k < octaves; ++k)
  {
    Array ridge_mask_current(z.shape);

    if (k > 0) apply_deposition = false;

    hmap::gpu::hydraulic_procedural(z,
                                    m * kp_global,
                                    a * c_erosion,
                                    seed,
                                    erosion_profile,
                                    erosion_profile_parameter,
                                    angle_shift,
                                    phase_smoothing,
                                    talus_ref,
                                    gradient_scaling_ratio,
                                    gradient_power,
                                    exclude_ridges,
                                    apply_deposition,
                                    deposition_strength,
                                    enable_default_noise,
                                    noise_amp,
                                    p_kp_multiplier,
                                    p_angle_shift,
                                    p_noise_x,
                                    p_noise_y,
                                    &ridge_mask_current,
                                    bbox);

    ridge_mask += a * ridge_mask_current;
    a_sum += a;

    a *= persistence;
    m *= lacunarity;
  }

  if (p_ridge_mask) *p_ridge_mask = 0.5f * (ridge_mask / a_sum) + 0.5f;
}

void hydraulic_procedural_fbm(Array         &z,
                              float          kp_global,
                              float          c_erosion,
                              std::uint32_t  seed,
                              const Array   *p_mask,
                              ErosionProfile erosion_profile,
                              int            octaves,
                              float          persistence,
                              float          lacunarity,
                              float          erosion_profile_parameter,
                              float          angle_shift,
                              float          phase_smoothing,
                              float          talus_ref,
                              float          gradient_scaling_ratio,
                              float          gradient_power,
                              bool           exclude_ridges,
                              bool           apply_deposition,
                              float          deposition_strength,
                              bool           enable_default_noise,
                              float          noise_amp,
                              const Array   *p_kp_multiplier,
                              const Array   *p_angle_shift,
                              const Array   *p_noise_x,
                              const Array   *p_noise_y,
                              Array         *p_ridge_mask,
                              glm::vec4      bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_kp_multiplier && !validate_same_shape(z, *p_kp_multiplier)) return;
  if (p_angle_shift && !validate_same_shape(z, *p_angle_shift)) return;
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return;
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return;

  if (!p_mask)
  {
    hydraulic_procedural_fbm(z,
                             kp_global,
                             c_erosion,
                             seed,
                             erosion_profile,
                             octaves,
                             persistence,
                             lacunarity,
                             erosion_profile_parameter,
                             angle_shift,
                             phase_smoothing,
                             talus_ref,
                             gradient_scaling_ratio,
                             gradient_power,
                             exclude_ridges,
                             apply_deposition,
                             deposition_strength,
                             enable_default_noise,
                             noise_amp,
                             p_kp_multiplier,
                             p_angle_shift,
                             p_noise_x,
                             p_noise_y,
                             p_ridge_mask,
                             bbox);
  }
  else
  {
    Array z_f = z;
    hydraulic_procedural_fbm(z_f,
                             kp_global,
                             c_erosion,
                             seed,
                             erosion_profile,
                             octaves,
                             persistence,
                             lacunarity,
                             erosion_profile_parameter,
                             angle_shift,
                             phase_smoothing,
                             talus_ref,
                             gradient_scaling_ratio,
                             gradient_power,
                             exclude_ridges,
                             apply_deposition,
                             deposition_strength,
                             enable_default_noise,
                             noise_amp,
                             p_kp_multiplier,
                             p_angle_shift,
                             p_noise_x,
                             p_noise_y,
                             p_ridge_mask,
                             bbox);
    z = lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap::gpu
