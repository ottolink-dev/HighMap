/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"

#include "highmap/functions.hpp"
#include "highmap/primitives/coherent_noise.hpp"

#include "NoiseLib/include/controlfunction.h"
#include "NoiseLib/include/math2d.h"
#include "NoiseLib/include/noise.h"

#pragma GCC diagnostic pop

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/internal/dendry_array_control_function.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

Array dendry(glm::ivec2    shape,
             glm::vec2     kw,
             std::uint32_t seed,
             Array        &control_array,
             float         eps,
             int           resolution,
             float         displacement,
             int           primitives_resolution_steps,
             float         slope_power,
             float         noise_amplitude_proportion,
             bool          add_control_function,
             float         control_function_overlap,
             const Array  *p_noise_x,
             const Array  *p_noise_y,
             glm::vec4     bbox,
             int           subsampling)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_non_empty(control_array)) return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  int nbuffer = (int)(control_function_overlap * control_array.shape.x);

  glm::ivec4 buffers = {nbuffer, nbuffer, nbuffer, nbuffer};
  Array      control_array_buffered = generate_buffered_array(control_array,
                                                         buffers);

  std::unique_ptr<ArrayControlFunction> control_function(
      std::make_unique<ArrayControlFunction>(control_array_buffered));

  const Point2D noise_top_left(0.f, 0.f);
  const Point2D noise_bottom_right(kw.x, kw.x);

  const Point2D control_function_top_left(0.5f * control_function_overlap,
                                          0.5f * control_function_overlap);
  const Point2D control_function_bottom_right(
      1.f - 0.5f * control_function_overlap,
      1.f - 0.5f * control_function_overlap);

  const Noise<ArrayControlFunction> noise(std::move(control_function),
                                          noise_top_left,
                                          noise_bottom_right,
                                          control_function_top_left,
                                          control_function_bottom_right,
                                          seed,
                                          eps,
                                          resolution,
                                          displacement,
                                          primitives_resolution_steps,
                                          slope_power,
                                          noise_amplitude_proportion,
                                          add_control_function,
                                          false,
                                          false,
                                          false,
                                          false);

  fill_array_using_xy_function(
      array,
      bbox,
      nullptr,
      p_noise_x,
      p_noise_y,
      nullptr,
      [&noise, &kw](float x, float y, float)
      { return noise.evaluateTerrain(kw.x * x, kw.y * y); },
      subsampling);

  return array;
}

Array dendry(glm::ivec2     shape,
             glm::vec2      kw,
             std::uint32_t  seed,
             NoiseFunction &noise_function,
             float          noise_function_offset,
             float          noise_function_scaling,
             float          eps,
             int            resolution,
             float          displacement,
             int            primitives_resolution_steps,
             float          slope_power,
             float          noise_amplitude_proportion,
             bool           add_control_function,
             float /* control_function_overlap */,
             const Array *p_noise_x,
             const Array *p_noise_y,
             glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  std::unique_ptr<XyControlFunction> control_function(
      std::make_unique<XyControlFunction>(noise_function,
                                          noise_function_offset,
                                          noise_function_scaling));

  const Point2D noise_top_left(0.f, 0.f);
  const Point2D noise_bottom_right(kw.x, kw.x);
  const Point2D control_function_top_left(0.f, 0.f);
  const Point2D control_function_bottom_right(1.f, 1.f);

  const Noise<XyControlFunction> noise(std::move(control_function),
                                       noise_top_left,
                                       noise_bottom_right,
                                       control_function_top_left,
                                       control_function_bottom_right,
                                       seed,
                                       eps,
                                       resolution,
                                       displacement,
                                       primitives_resolution_steps,
                                       slope_power,
                                       noise_amplitude_proportion,
                                       add_control_function,
                                       false,
                                       false,
                                       false,
                                       false);

  fill_array_using_xy_function(
      array,
      bbox,
      nullptr,
      p_noise_x,
      p_noise_y,
      nullptr,
      [&noise, &kw](float x, float y, float)
      { return noise.evaluateTerrain(kw.x * x, kw.y * y); });

  return array;
}

} // namespace hmap
