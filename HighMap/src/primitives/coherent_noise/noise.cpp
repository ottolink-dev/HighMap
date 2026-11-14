/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <memory>
#include <utility>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/coherent_noise.hpp"

namespace hmap
{

Array noise(NoiseType     noise_type,
            glm::ivec2    shape,
            glm::vec2     kw,
            std::uint32_t seed,
            const Array  *p_noise_x,
            const Array  *p_noise_y,
            glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  with_noise_function(noise_type,
                      kw,
                      seed,
                      [&](const auto &base_fn)
                      {
                        fill_array_using_xy_function(array,
                                                     bbox,
                                                     nullptr,
                                                     p_noise_x,
                                                     p_noise_y,
                                                     nullptr,
                                                     base_fn);
                      });
  return array;
}

Array noise_fbm(NoiseType     noise_type,
                glm::ivec2    shape,
                glm::vec2     kw,
                std::uint32_t seed,
                int           octaves,
                float         weight,
                float         persistence,
                float         lacunarity,
                const Array  *p_ctrl_param,
                const Array  *p_noise_x,
                const Array  *p_noise_y,
                glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  with_noise_function(noise_type,
                      kw,
                      seed,
                      [&](const auto &base_fn)
                      {
                        hmap::FbmFunction f = hmap::FbmFunction(base_fn,
                                                                octaves,
                                                                weight,
                                                                persistence,
                                                                lacunarity);
                        fill_array_using_xy_function(array,
                                                     bbox,
                                                     p_ctrl_param,
                                                     p_noise_x,
                                                     p_noise_y,
                                                     nullptr,
                                                     f);
                      });
  return array;
}

Array noise_iq(NoiseType     noise_type,
               glm::ivec2    shape,
               glm::vec2     kw,
               std::uint32_t seed,
               int           octaves,
               float         weight,
               float         persistence,
               float         lacunarity,
               float         gradient_scale,
               const Array  *p_ctrl_param,
               const Array  *p_noise_x,
               const Array  *p_noise_y,
               glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  with_noise_function(noise_type,
                      kw,
                      seed,
                      [&](const auto &base_fn)
                      {
                        hmap::FbmIqFunction f = hmap::FbmIqFunction(
                            base_fn,
                            octaves,
                            weight,
                            persistence,
                            lacunarity,
                            gradient_scale);
                        fill_array_using_xy_function(array,
                                                     bbox,
                                                     p_ctrl_param,
                                                     p_noise_x,
                                                     p_noise_y,
                                                     nullptr,
                                                     f);
                      });
  return array;
}

Array noise_jordan(NoiseType     noise_type,
                   glm::ivec2    shape,
                   glm::vec2     kw,
                   std::uint32_t seed,
                   int           octaves,
                   float         weight,
                   float         persistence,
                   float         lacunarity,
                   float         warp0,
                   float         damp0,
                   float         warp_scale,
                   float         damp_scale,
                   const Array  *p_ctrl_param,
                   const Array  *p_noise_x,
                   const Array  *p_noise_y,
                   glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  with_noise_function(noise_type,
                      kw,
                      seed,
                      [&](const auto &base_fn)
                      {
                        hmap::FbmJordanFunction f = hmap::FbmJordanFunction(
                            base_fn,
                            octaves,
                            weight,
                            persistence,
                            lacunarity,
                            warp0,
                            damp0,
                            warp_scale,
                            damp_scale);
                        fill_array_using_xy_function(array,
                                                     bbox,
                                                     p_ctrl_param,
                                                     p_noise_x,
                                                     p_noise_y,
                                                     nullptr,
                                                     f);
                      });
  return array;
}

Array noise_parberry(glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     int           octaves,
                     float         weight,
                     float         persistence,
                     float         lacunarity,
                     float         mu,
                     const Array  *p_ctrl_param,
                     const Array  *p_noise_x,
                     const Array  *p_noise_y,
                     glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  std::unique_ptr<NoiseFunction> p = std::unique_ptr<NoiseFunction>(
      new ParberryFunction(kw, seed, mu));

  hmap::FbmFunction f = hmap::FbmFunction(ParberryFunction(kw, seed, mu),
                                          octaves,
                                          weight,
                                          persistence,
                                          lacunarity);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f);
  return array;
}

Array noise_pingpong(NoiseType     noise_type,
                     glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     int           octaves,
                     float         weight,
                     float         persistence,
                     float         lacunarity,
                     const Array  *p_ctrl_param,
                     const Array  *p_noise_x,
                     const Array  *p_noise_y,
                     glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  with_noise_function(noise_type,
                      kw,
                      seed,
                      [&](const auto &base_fn)
                      {
                        hmap::FbmPingpongFunction f = hmap::FbmPingpongFunction(
                            base_fn,
                            octaves,
                            weight,
                            persistence,
                            lacunarity);
                        fill_array_using_xy_function(array,
                                                     bbox,
                                                     p_ctrl_param,
                                                     p_noise_x,
                                                     p_noise_y,
                                                     nullptr,
                                                     f);
                      });
  return array;
}

Array noise_ridged(NoiseType     noise_type,
                   glm::ivec2    shape,
                   glm::vec2     kw,
                   std::uint32_t seed,
                   int           octaves,
                   float         weight,
                   float         persistence,
                   float         lacunarity,
                   float         k_smoothing,
                   const Array  *p_ctrl_param,
                   const Array  *p_noise_x,
                   const Array  *p_noise_y,
                   glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  with_noise_function(noise_type,
                      kw,
                      seed,
                      [&](const auto &base_fn)
                      {
                        hmap::FbmRidgedFunction f = hmap::FbmRidgedFunction(
                            base_fn,
                            octaves,
                            weight,
                            persistence,
                            lacunarity,
                            k_smoothing);
                        fill_array_using_xy_function(array,
                                                     bbox,
                                                     p_ctrl_param,
                                                     p_noise_x,
                                                     p_noise_y,
                                                     nullptr,
                                                     f);
                      });
  return array;
}

Array noise_swiss(NoiseType     noise_type,
                  glm::ivec2    shape,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  int           octaves,
                  float         weight,
                  float         persistence,
                  float         lacunarity,
                  float         warp_scale,
                  const Array  *p_ctrl_param,
                  const Array  *p_noise_x,
                  const Array  *p_noise_y,
                  glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  with_noise_function(noise_type,
                      kw,
                      seed,
                      [&](const auto &base_fn)
                      {
                        hmap::FbmSwissFunction f = hmap::FbmSwissFunction(
                            base_fn,
                            octaves,
                            weight,
                            persistence,
                            lacunarity,
                            warp_scale);
                        fill_array_using_xy_function(array,
                                                     bbox,
                                                     p_ctrl_param,
                                                     p_noise_x,
                                                     p_noise_y,
                                                     nullptr,
                                                     f);
                      });
  return array;
}

} // namespace hmap
