/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

void phase_averaging(Array &field_real, Array &field_imag, int ir)
{
  if (!validate_non_empty(field_real) ||
      !validate_same_shape(field_real, field_imag))
  {
    return;
  }

  const glm::ivec2 shape = field_real.shape;

  auto run = clwrapper::Run("phase_averaging");

  // inputs
  run.bind_imagef("fr", field_real.vector, shape.x, shape.y);
  run.bind_imagef("fi", field_imag.vector, shape.x, shape.y);

  // outputs
  run.bind_imagef("fr_out", field_real.vector, shape.x, shape.y, true);
  run.bind_imagef("fi_out", field_imag.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x, shape.y, ir);

  run.execute({shape.x, shape.y});

  // update flux (from GPU to CPU)
  run.read_imagef("fr_out");
  run.read_imagef("fi_out");
}

Array phase_field(const Array     &array,
                  const glm::vec2 &kw,
                  std::uint32_t    seed,
                  float            kp,
                  bool             rotate90,
                  float            normalization,
                  const glm::vec2 &jitter,
                  int              angle_filter_ir,
                  const Array     *p_ctrl_param,
                  const Array     *p_noise_x,
                  const Array     *p_noise_y,
                  Array           *p_modulus,
                  Array           *p_angle_jump_mask,
                  glm::vec4        bbox)
{
  if (!validate_non_empty(array)) return Array();
  if (p_ctrl_param && !validate_same_shape(array, *p_ctrl_param))
    return Array(array.shape);
  if (p_noise_x && !validate_same_shape(array, *p_noise_x))
    return Array(array.shape);
  if (p_noise_y && !validate_same_shape(array, *p_noise_y))
    return Array(array.shape);

  const glm::ivec2 shape = array.shape;

  // --- compute local angle

  float phi = rotate90 ? M_PI : 0.5 * M_PI;
  Array dx = gradient_x(array);
  Array dy = gradient_y(array);
  phase_averaging(dx, dy, angle_filter_ir);
  Array angle = atan2(dy, dx) + phi;

  if (p_angle_jump_mask)
    *p_angle_jump_mask = talus_jump_mask(angle, M_PI, 0.1f * M_PI);

  // --- compute phase

  Array phase(shape);

  auto run = clwrapper::Run("phase_field");

  run.bind_buffer<float>("angle", angle.vector);
  run.bind_buffer<float>("phase", phase.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "p_noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "p_noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "p_modulus", p_modulus);

  run.bind_arguments(shape.x,
                     shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     jitter,
                     normalization,
                     kp,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_modulus ? 1 : 0,
                     bbox);

  run.write_buffer("angle");
  run.write_buffer("phase");

  run.execute({shape.x, shape.y});

  run.read_buffer("phase");
  if (p_modulus) run.read_buffer("p_modulus");

  return phase;
}

Array phase_field(const Array     &array,
                  std::uint32_t    seed,
                  float            kp_global,
                  bool             rotate90,
                  float            normalization,
                  const glm::vec2 &jitter,
                  int              angle_filter_ir,
                  const Array     *p_ctrl_param,
                  const Array     *p_noise_x,
                  const Array     *p_noise_y,
                  Array           *p_modulus,
                  Array           *p_angle_jump_mask,
                  glm::vec4        bbox)
{
  if (!validate_non_empty(array)) return Array();
  if (p_ctrl_param && !validate_same_shape(array, *p_ctrl_param))
    return Array(array.shape);
  if (p_noise_x && !validate_same_shape(array, *p_noise_x))
    return Array(array.shape);
  if (p_noise_y && !validate_same_shape(array, *p_noise_y))
    return Array(array.shape);

  float           kp = 1.f;
  const glm::vec2 kw = {kp_global, kp_global};

  return phase_field(array,
                     kw,
                     seed,
                     kp,
                     rotate90,
                     normalization,
                     jitter,
                     angle_filter_ir,
                     p_ctrl_param,
                     p_noise_x,
                     p_noise_y,
                     p_modulus,
                     p_angle_jump_mask,
                     bbox);
}

Array phase_field_angle(const Array     &angle,
                        const glm::vec2 &kw,
                        std::uint32_t    seed,
                        float            kp,
                        float            normalization,
                        const glm::vec2 &jitter,
                        const Array     *p_ctrl_param,
                        const Array     *p_noise_x,
                        const Array     *p_noise_y,
                        Array           *p_modulus,
                        Array           *p_angle_jump_mask,
                        glm::vec4        bbox)
{
  if (!validate_non_empty(angle)) return Array();
  if (p_ctrl_param && !validate_same_shape(angle, *p_ctrl_param))
    return Array(angle.shape);
  if (p_noise_x && !validate_same_shape(angle, *p_noise_x))
    return Array(angle.shape);
  if (p_noise_y && !validate_same_shape(angle, *p_noise_y))
    return Array(angle.shape);

  const glm::ivec2 shape = angle.shape;

  if (p_angle_jump_mask)
    *p_angle_jump_mask = talus_jump_mask(angle, M_PI, 0.1f * M_PI);

  // --- compute phase

  Array phase(shape);

  auto run = clwrapper::Run("phase_field");

  run.bind_buffer<float>("angle", angle.vector);
  run.bind_buffer<float>("phase", phase.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "p_noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "p_noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "p_modulus", p_modulus);

  run.bind_arguments(shape.x,
                     shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     jitter,
                     normalization,
                     kp,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_modulus ? 1 : 0,
                     bbox);

  run.write_buffer("angle");
  run.write_buffer("phase");

  run.execute({shape.x, shape.y});

  run.read_buffer("phase");
  if (p_modulus) run.read_buffer("p_modulus");

  return phase;
}

} // namespace hmap::gpu
