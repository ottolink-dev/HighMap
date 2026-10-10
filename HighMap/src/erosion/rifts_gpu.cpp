/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

void rifts(Array           &z,
           const glm::vec2 &kw,
           float            angle,
           float            amplitude,
           std::uint32_t    seed,
           float            elevation_noise_shift,
           float            k_smooth_bottom,
           float            k_smooth_top,
           float            radial_spread_amp,
           float            elevation_noise_amp,
           float            clamp_vmin,
           float            remap_vmin,
           bool             apply_mask,
           bool             reverse_mask,
           float            mask_gamma,
           const Array     *p_noise_x,
           const Array     *p_noise_y,
           const Array     *p_mask,
           const glm::vec2 &center,
           const glm::vec4 &bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return;
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  auto run = clwrapper::Run("rifts");

  run.bind_buffer<float>("z", z.vector);

  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "mask", p_mask);

  run.bind_arguments(z.shape.x,
                     z.shape.y,
                     kw,
                     angle,
                     amplitude,
                     seed,
                     elevation_noise_shift,
                     k_smooth_bottom,
                     k_smooth_top,
                     radial_spread_amp,
                     elevation_noise_amp,
                     clamp_vmin,
                     remap_vmin,
                     apply_mask ? 1 : 0,
                     reverse_mask ? 1 : 0,
                     mask_gamma,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_mask ? 1 : 0,
                     center,
                     bbox);

  run.write_buffer("z");
  run.execute({z.shape.x, z.shape.y});
  run.read_buffer("z");
}

} // namespace hmap::gpu
