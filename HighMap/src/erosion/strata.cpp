/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/erosion/strata_erosion.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/coherent_noise.hpp"

namespace hmap::gpu
{

void strata(Array           &z,
            float            angle,
            float            slope,
            float            gamma,
            std::uint32_t    seed,
            bool             linear_gamma,
            float            kz,
            int              octaves,
            float            lacunarity,
            float            gamma_noise_ratio,
            float            noise_amp,
            const glm::vec2 &noise_kw,
            bool             enable_ridge_noise,
            const glm::vec2 &ridge_noise_kw,
            float            ridge_angle_shift,
            float            ridge_noise_amp,
            float            ridge_clamp_vmin,
            float            ridge_remap_vmin,
            bool             apply_elevation_mask,
            bool             apply_ridge_mask,
            float            mask_gamma,
            const Array     *p_mask,
            const glm::vec4 &bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  auto run = clwrapper::Run("strata");

  run.bind_buffer<float>("z", z.vector);

  helper_bind_optional_buffer(run, "mask", p_mask);

  run.bind_arguments(z.shape.x,
                     z.shape.y,
                     angle,
                     slope,
                     gamma,
                     seed,
                     linear_gamma ? 1 : 0,
                     kz,
                     octaves,
                     lacunarity,
                     gamma_noise_ratio,
                     noise_amp,
                     noise_kw,
                     enable_ridge_noise ? 1 : 0,
                     ridge_noise_kw,
                     ridge_angle_shift,
                     ridge_noise_amp,
                     ridge_clamp_vmin,
                     ridge_remap_vmin,
                     apply_elevation_mask ? 1 : 0,
                     apply_ridge_mask ? 1 : 0,
                     mask_gamma,
                     p_mask ? 1 : 0,
                     bbox);

  run.write_buffer("z");
  run.execute({z.shape.x, z.shape.y});
  run.read_buffer("z");
}

void strata_cells(Array        &z,
                  glm::vec2     kw,
                  float         amp,
                  std::uint32_t seed,
                  float         gamma,
                  float         gamma_lateral,
                  float         angle,
                  float         noise_amp,
                  bool          absolute_displacement,
                  float         occurence_probability,
                  const Array  *p_noise_x,
                  const Array  *p_noise_y,
                  glm::vec4     bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return;
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return;

  const Array *p_mask = nullptr;

  //

  auto run = clwrapper::Run("strata_cells");

  run.bind_buffer<float>("z", z.vector);

  helper_bind_optional_buffer(run, "mask", p_mask);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(z.shape.x,
                     z.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     gamma,
                     gamma_lateral,
                     angle,
                     amp,
                     noise_amp,
                     absolute_displacement ? 1 : 0,
                     occurence_probability,
                     p_mask ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_buffer("z");
  run.execute({z.shape.x, z.shape.y});
  run.read_buffer("z");
}

void strata_cells(Array        &z,
                  glm::vec2     kw,
                  float         amp,
                  std::uint32_t seed,
                  const Array  *p_mask,
                  float         gamma,
                  float         gamma_lateral,
                  float         angle,
                  float         noise_amp,
                  bool          absolute_displacement,
                  float         occurence_probability,
                  const Array  *p_noise_x,
                  const Array  *p_noise_y,
                  glm::vec4     bbox)
{
  apply_with_mask(z,
                  p_mask,
                  [&](Array &a)
                  {
                    strata_cells(a,
                                 kw,
                                 amp,
                                 seed,
                                 gamma,
                                 gamma_lateral,
                                 angle,
                                 noise_amp,
                                 absolute_displacement,
                                 occurence_probability,
                                 p_noise_x,
                                 p_noise_y,
                                 bbox);
                  });
}

void strata_cells_fbm(Array        &z,
                      glm::vec2     kw,
                      float         amp,
                      std::uint32_t seed,
                      float         gamma,
                      float         gamma_lateral,
                      float         angle,
                      bool          enable_default_noise,
                      float         default_noise_amp,
                      bool          absolute_displacement,
                      float         occurence_probability,
                      int           octaves,
                      float         persistence,
                      float         lacunarity,
                      const Array  *p_noise_x,
                      const Array  *p_noise_y,
                      glm::vec4     bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return;
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return;

  // --- Default noise

  Array dx, dy;

  if (enable_default_noise)
  {
    const glm::vec2 kw_noise = kw;

    dx = default_noise_amp * gpu::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                            z.shape,
                                            kw_noise,
                                            ++seed,
                                            /* octaves */ 8,
                                            /* weight*/ 0.f,
                                            /* persistence */ 0.5f,
                                            /* lacunarity */ 2.f,
                                            nullptr,
                                            nullptr,
                                            nullptr,
                                            bbox);

    dy = default_noise_amp * gpu::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                            z.shape,
                                            kw_noise,
                                            ++seed,
                                            /* octaves */ 8,
                                            /* weight*/ 0.f,
                                            /* persistence */ 0.5f,
                                            /* lacunarity */ 2.f,
                                            nullptr,
                                            nullptr,
                                            nullptr,
                                            bbox);

    p_noise_x = &dx;
    p_noise_y = &dy;
  }

  // --- Fbm layers

  float a = 1.f;
  float m = 1.f;

  for (int k = 0; k < octaves; ++k)
  {
    strata_cells(z,
                 m * kw,
                 a * amp,
                 seed,
                 gamma,
                 gamma_lateral,
                 angle,
                 /* noise_amp */ 1.f,
                 absolute_displacement,
                 occurence_probability,
                 p_noise_x,
                 p_noise_y,
                 bbox);

    a *= persistence;
    m *= lacunarity;
  }
}

void strata_cells_fbm(Array        &z,
                      glm::vec2     kw,
                      float         amp,
                      std::uint32_t seed,
                      const Array  *p_mask,
                      float         gamma,
                      float         gamma_lateral,
                      float         angle,
                      bool          enable_default_noise,
                      float         default_noise_amp,
                      bool          absolute_displacement,
                      float         occurence_probability,
                      int           octaves,
                      float         persistence,
                      float         lacunarity,
                      const Array  *p_noise_x,
                      const Array  *p_noise_y,
                      glm::vec4     bbox)
{
  apply_with_mask(z,
                  p_mask,
                  [&](Array &a)
                  {
                    strata_cells_fbm(a,
                                     kw,
                                     amp,
                                     seed,
                                     gamma,
                                     gamma_lateral,
                                     angle,
                                     enable_default_noise,
                                     default_noise_amp,
                                     absolute_displacement,
                                     occurence_probability,
                                     octaves,
                                     persistence,
                                     lacunarity,
                                     p_noise_x,
                                     p_noise_y,
                                     bbox);
                  });
}

void strata_terrace(Array        &z,
                    float         gamma,
                    std::uint32_t seed,
                    float         kz,
                    bool          linear_gamma,
                    float         gamma_noise_ratio,
                    float         slope,
                    float         angle,
                    const Array  *p_noise,
                    glm::vec4     bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_noise && !validate_same_shape(z, *p_noise)) return;

  auto run = clwrapper::Run("strata_terrace");

  run.bind_buffer<float>("z", z.vector);
  helper_bind_optional_buffer(run, "noise", p_noise);

  run.bind_arguments(z.shape.x,
                     z.shape.y,
                     gamma,
                     seed,
                     linear_gamma ? 1 : 0,
                     kz,
                     gamma_noise_ratio,
                     slope,
                     angle,
                     p_noise ? 1 : 0,
                     bbox);

  run.write_buffer("z");
  run.execute({z.shape.x, z.shape.y});
  run.read_buffer("z");
}

void strata_terrace(Array        &z,
                    float         gamma,
                    std::uint32_t seed,
                    const Array  *p_mask,
                    float         kz,
                    bool          linear_gamma,
                    float         gamma_noise_ratio,
                    float         slope,
                    float         angle,
                    const Array  *p_noise,
                    glm::vec4     bbox)
{
  apply_with_mask(z,
                  p_mask,
                  [&](Array &a)
                  {
                    strata_terrace(a,
                                   gamma,
                                   seed,
                                   kz,
                                   linear_gamma,
                                   gamma_noise_ratio,
                                   slope,
                                   angle,
                                   p_noise,
                                   bbox);
                  });
}

} // namespace hmap::gpu
