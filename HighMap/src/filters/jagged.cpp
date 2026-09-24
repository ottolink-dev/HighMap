/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

Array jagged(const Array  &array,
             glm::vec2     kw,
             float         amp,
             std::uint32_t seed,
             glm::vec2     jitter,
             float         gamma,
             float         angle,
             const Array  *p_mask,
             const Array  *p_noise_x,
             const Array  *p_noise_y,
             glm::vec4     bbox)
{
  if (!validate_non_empty(array)) return Array();
  if (p_mask && !validate_same_shape(array, *p_mask)) return Array(array.shape);
  if (p_noise_x && !validate_same_shape(array, *p_noise_x))
    return Array(array.shape);
  if (p_noise_y && !validate_same_shape(array, *p_noise_y))
    return Array(array.shape);

  const glm::ivec2 shape = array.shape;
  Array            out(shape);

  auto run = clwrapper::Run("jagged");

  run.bind_imagef("in", array.vector, shape.x, shape.y);
  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  helper_bind_optional_buffer(run, "mask", p_mask);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(shape.x,
                     shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     jitter,
                     amp,
                     gamma,
                     angle,
                     p_mask ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_imagef("in");
  run.execute({shape.x, shape.y});
  run.read_imagef("out");

  return out;
}

Array jagged(const Array  &array,
             float         kw,
             float         amp,
             std::uint32_t seed,
             glm::vec2     jitter,
             float         gamma,
             float         angle,
             const Array  *p_mask,
             const Array  *p_noise_x,
             const Array  *p_noise_y,
             glm::vec4     bbox)
{
  return jagged(array,
                {kw, kw},
                amp,
                seed,
                jitter,
                gamma,
                angle,
                p_mask,
                p_noise_x,
                p_noise_y,
                bbox);
}

Array jagged_fbm(const Array  &array,
                 glm::vec2     kw,
                 float         amp,
                 std::uint32_t seed,
                 int           octaves,
                 float         persistence,
                 float         lacunarity,
                 bool          switch_kx_ky,
                 glm::vec2     jitter,
                 float         gamma,
                 float         angle,
                 const Array  *p_mask,
                 const Array  *p_noise_x,
                 const Array  *p_noise_y,
                 glm::vec4     bbox)
{
  if (!validate_non_empty(array)) return Array();
  if (p_mask && !validate_same_shape(array, *p_mask)) return Array(array.shape);
  if (p_noise_x && !validate_same_shape(array, *p_noise_x))
    return Array(array.shape);
  if (p_noise_y && !validate_same_shape(array, *p_noise_y))
    return Array(array.shape);

  Array         out = array;
  glm::vec2     current_kw = kw;
  float         current_amp = amp;
  std::uint32_t current_seed = seed;

  for (int k = 0; k < octaves; ++k)
  {
    out = jagged(out,
                 current_kw,
                 current_amp,
                 current_seed,
                 jitter,
                 gamma,
                 angle,
                 p_mask,
                 p_noise_x,
                 p_noise_y,
                 bbox);

    current_kw *= lacunarity;
    if (switch_kx_ky) std::swap(current_kw.x, current_kw.y);

    current_amp *= persistence;
    current_seed++;
  }

  return out;
}

Array jagged_fbm(const Array  &array,
                 float         kw,
                 float         amp,
                 std::uint32_t seed,
                 int           octaves,
                 float         persistence,
                 float         lacunarity,
                 bool          switch_kx_ky,
                 glm::vec2     jitter,
                 float         gamma,
                 float         angle,
                 const Array  *p_mask,
                 const Array  *p_noise_x,
                 const Array  *p_noise_y,
                 glm::vec4     bbox)
{
  return jagged_fbm(array,
                    {kw, kw},
                    amp,
                    seed,
                    octaves,
                    persistence,
                    lacunarity,
                    switch_kx_ky,
                    jitter,
                    gamma,
                    angle,
                    p_mask,
                    p_noise_x,
                    p_noise_y,
                    bbox);
}

} // namespace hmap::gpu
