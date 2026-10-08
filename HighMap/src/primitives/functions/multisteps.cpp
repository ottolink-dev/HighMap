/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>
#include <functional>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap
{

Array multisteps(glm::ivec2       shape,
                 float            angle,
                 float            r,
                 int              nsteps,
                 float            elevation_exponent,
                 float            shape_gain,
                 float            scale,
                 float            outer_slope,
                 const Array     *p_ctrl_param,
                 const Array     *p_noise_x,
                 const Array     *p_noise_y,
                 const glm::vec2 &center,
                 const glm::vec4 &bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  float alpha = angle / 180.f * M_PI;
  float ca = std::cos(alpha);
  float sa = std::sin(alpha);
  float d0 = r == 1.f ? 1.f / float(nsteps)
                      : (1.f - r) / (1.f - std::pow(r, nsteps));

  auto lambda = [=](float x, float y, float ctrl)
  {
    float t = (x - center.x + 0.5f * scale) * ca +
              (y - center.y + 0.5f * scale) * sa;
    t /= scale;

    // out-of-range
    if (t < 0.f)
      return 1.f - outer_slope * t;
    else if (t > 1.f)
      return outer_slope * (1.f - t);

    // find correspond step
    int   n = std::floor(std::log((d0 + (r - 1.f) * t) / d0) / std::log(r));
    float ts = d0 * (1.f - std::pow(r, n)) / (1.f - r);
    float te = d0 * (1.f - std::pow(r, n + 1)) / (1.f - r);

    float zs = 1.f - std::pow(float(n) / float(nsteps), elevation_exponent);
    float ze = 1.f - std::pow(float(n + 1) / float(nsteps), elevation_exponent);

    float s = (t - ts) / (te - ts); // in [0, 1]
    s = std::pow(s, shape_gain);

    float v = lerp(zs, ze, s); // steps
    float v_linear = 1.f - t;  // linear transition

    // control parameters control balance between steps and linear
    // transition (by default ctrl = 1)
    return lerp(v, v_linear, 1.f - ctrl);
  };

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               lambda);

  return array;
}

} // namespace hmap

namespace hmap::gpu
{

Array multisteps(glm::ivec2       shape,
                 float            angle,
                 std::uint32_t    seed,
                 glm::vec2        kw,
                 float            noise_amp,
                 float            noise_rugosity,
                 bool             noise_inflate,
                 float            r,
                 int              nsteps,
                 float            elevation_exponent,
                 float            shape_gain,
                 float            scale,
                 float            outer_slope,
                 const Array     *p_ctrl_param,
                 const glm::vec2 &center,
                 const glm::vec4 &bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);

  // built-in noise
  glm::vec2         jitter(1.f, 1.f);
  VoronoiReturnType return_type = VoronoiReturnType::EDGE_DISTANCE_SQUARED;

  Array noise = 2.f * gpu::voronoi_fbm(shape,
                                       kw,
                                       seed,
                                       jitter,
                                       /* bias */ 0.f,
                                       /* k_smoothing */ 0.f,
                                       0.f,
                                       return_type,
                                       /* octaves */ 8,
                                       noise_rugosity,
                                       /* persistence */ 0.5f,
                                       /* lacunarity */ 2.f,
                                       nullptr,
                                       nullptr,
                                       nullptr,
                                       bbox); // in [0, 1]

  float sign = noise_inflate ? -1.f : 1.f;
  noise = noise_amp * sign * (2.f * noise - 1.f);

  float alpha = angle / 180.f * M_PI;
  Array dx = std::cos(alpha) * noise;
  Array dy = std::sin(alpha) * noise;

  return multisteps(shape,
                    angle,
                    r,
                    nsteps,
                    elevation_exponent,
                    shape_gain,
                    scale,
                    outer_slope,
                    p_ctrl_param,
                    &dx,
                    &dy,
                    center,
                    bbox);
}

} // namespace hmap::gpu
