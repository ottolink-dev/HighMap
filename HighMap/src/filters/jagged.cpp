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
             std::uint32_t seed,
             glm::vec2     jitter,
             float         gamma,
             float         shape_gamma,
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
                     gamma,
                     shape_gamma,
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
             std::uint32_t seed,
             glm::vec2     jitter,
             float         gamma,
             float         shape_gamma,
             const Array  *p_mask,
             const Array  *p_noise_x,
             const Array  *p_noise_y,
             glm::vec4     bbox)
{
  return jagged(array,
                {kw, kw},
                seed,
                jitter,
                gamma,
                shape_gamma,
                p_mask,
                p_noise_x,
                p_noise_y,
                bbox);
}

} // namespace hmap::gpu
