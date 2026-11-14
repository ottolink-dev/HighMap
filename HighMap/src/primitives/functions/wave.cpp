/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap
{

Array wave_dune(glm::ivec2   shape,
                float        kw,
                float        angle,
                float        xtop,
                float        xbottom,
                float        phase_shift,
                const Array *p_noise_x,
                const Array *p_noise_y,
                glm::vec2    center,
                glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array                  array = Array(shape);
  hmap::WaveDuneFunction f = hmap::WaveDuneFunction({kw, kw},
                                                    angle,
                                                    xtop,
                                                    xbottom,
                                                    phase_shift,
                                                    center);

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f);
  return array;
}

Array wave_sine(glm::ivec2   shape,
                float        kw,
                float        angle,
                float        phase_shift,
                const Array *p_noise_x,
                const Array *p_noise_y,
                glm::vec2    center,
                glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array                  array = Array(shape);
  hmap::WaveSineFunction f = hmap::WaveSineFunction({kw, kw},
                                                    angle,
                                                    phase_shift,
                                                    center);

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f);
  return array;
}

Array wave_square(glm::ivec2   shape,
                  float        kw,
                  float        angle,
                  float        phase_shift,
                  const Array *p_noise_x,
                  const Array *p_noise_y,
                  glm::vec2    center,
                  glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array                    array = Array(shape);
  hmap::WaveSquareFunction f = hmap::WaveSquareFunction({kw, kw},
                                                        angle,
                                                        phase_shift,
                                                        center);

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f);
  return array;
}

Array wave_triangular(glm::ivec2   shape,
                      float        kw,
                      float        angle,
                      float        slant_ratio,
                      float        phase_shift,
                      const Array *p_noise_x,
                      const Array *p_noise_y,
                      glm::vec2    center,
                      glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array                        array = Array(shape);
  hmap::WaveTriangularFunction f = hmap::WaveTriangularFunction({kw, kw},
                                                                angle,
                                                                slant_ratio,
                                                                phase_shift,
                                                                center);

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f);
  return array;
}

} // namespace hmap
