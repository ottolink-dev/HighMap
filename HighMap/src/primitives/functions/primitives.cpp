/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap
{

Array biquad_pulse(glm::ivec2   shape,
                   float        gain,
                   const Array *p_ctrl_param,
                   const Array *p_noise_x,
                   const Array *p_noise_y,
                   glm::vec2    center,
                   glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array          array = Array(shape);
  BiquadFunction f = BiquadFunction(gain, center);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

Array biquad_pulse_x(glm::ivec2 shape, glm::vec4 bbox)
{
  if (!validate_shape(shape)) return Array();

  Array array = Array(shape);

  auto fct = [](float x, float, float)
  {
    x = 16.f * x * (1.f - x);
    return std::clamp(x, 0.f, 1.f);
  };

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               nullptr,
                               nullptr,
                               nullptr,
                               fct);
  return array;
}

Array biquad_pulse_y(glm::ivec2 shape, glm::vec4 bbox)
{
  if (!validate_shape(shape)) return Array();

  Array array = Array(shape);

  auto fct = [](float, float y, float)
  {
    y = 16.f * y * (1.f - y);
    return std::clamp(y, 0.f, 1.f);
  };

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               nullptr,
                               nullptr,
                               nullptr,
                               fct);
  return array;
}

Array bump(glm::ivec2   shape,
           float        gain,
           const Array *p_ctrl_param,
           const Array *p_noise_x,
           const Array *p_noise_y,
           glm::vec2    center,
           glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array        array = Array(shape);
  BumpFunction f = BumpFunction(gain, center);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

Array bump_lorentzian(glm::ivec2   shape,
                      float        width_factor,
                      float        radius,
                      const Array *p_ctrl_param,
                      const Array *p_noise_x,
                      const Array *p_noise_y,
                      glm::vec2    center,
                      glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  float radius_sq = radius * radius;

  auto lambda =
      [radius_sq, width_factor, center](float x, float y, float ctrl_param)
  {
    float dx = x - center.x;
    float dy = y - center.y;
    float r2 = (dx * dx + dy * dy) / radius_sq;

    float width_sq = width_factor * width_factor * ctrl_param * ctrl_param;
    float c = 1.f / (1.f + 1.f / width_sq); // normalization coeff

    if (r2 < 1.f)
      return (1.f / (1.f + r2 / width_sq) - c) / (1.f - c);
    else
      return 0.f;
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

Array constant(glm::ivec2 shape, float value)
{
  if (!validate_shape(shape)) return Array();

  Array array = Array(shape);
  for (auto &v : array.vector)
    v = value;
  return array;
}

Array cubic_pulse(glm::ivec2   shape,
                  const Array *p_noise_x,
                  const Array *p_noise_y,
                  glm::vec2    center,
                  glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  auto lambda = [center](float x, float y, float)
  {
    float dx = x - center.x;
    float dy = y - center.y;
    float r = std::hypot(dx, dy) / 0.5f;

    if (r < 1.f)
      return 1.f - r * r * (3.f - 2.f * r);
    else
      return 0.f;
  };

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               lambda);
  return array;
}

Array disk(glm::ivec2   shape,
           float        radius,
           float        slope,
           const Array *p_ctrl_param,
           const Array *p_noise_x,
           const Array *p_noise_y,
           glm::vec2    center,
           glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array        array = Array(shape);
  DiskFunction f = DiskFunction(radius, slope, center);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

Array gaussian_pulse(glm::ivec2   shape,
                     float        sigma,
                     const Array *p_ctrl_param,
                     const Array *p_noise_x,
                     const Array *p_noise_y,
                     glm::vec2    center,
                     glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array                 array = Array(shape);
  GaussianPulseFunction f = GaussianPulseFunction(sigma, center);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

Array paraboloid(glm::ivec2   shape,
                 float        angle,
                 float        a,
                 float        b,
                 float        v0,
                 bool         reverse_x,
                 bool         reverse_y,
                 const Array *p_noise_x,
                 const Array *p_noise_y,
                 glm::vec2    center,
                 glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  float ca = std::cos(-angle / 180.f * M_PI);
  float sa = std::sin(-angle / 180.f * M_PI);

  float inv_a2 = (reverse_x ? -1.f : 1.f) * 1.f / (a * a);
  float inv_b2 = (reverse_y ? -1.f : 1.f) * 1.f / (b * b);

  auto lambda =
      [&ca, &sa, &v0, &inv_a2, &inv_b2, &center](float x, float y, float)
  {
    float xr = ca * (x - center.x) - sa * (y - center.y);
    float yr = sa * (x - center.x) + ca * (y - center.y);

    return inv_a2 * xr * xr + inv_b2 * yr * yr + v0;
  };

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               lambda);
  return array;
}

Array quad_surface(glm::ivec2   shape,
                   float        c00,
                   float        c10,
                   float        c01,
                   float        c11,
                   const Array *p_ctrl_param,
                   const Array *p_noise_x,
                   const Array *p_noise_y,
                   glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array               array = Array(shape);
  QuadSurfaceFunction f = QuadSurfaceFunction(c00, c10, c01, c11, bbox);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

Array rectangle(glm::ivec2   shape,
                float        rx,
                float        ry,
                float        angle,
                float        slope,
                const Array *p_ctrl_param,
                const Array *p_noise_x,
                const Array *p_noise_y,
                glm::vec2    center,
                glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array             array = Array(shape);
  RectangleFunction f = RectangleFunction(rx, ry, angle, slope, center);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

Array slope(glm::ivec2   shape,
            float        angle,
            float        slope,
            const Array *p_ctrl_param,
            const Array *p_noise_x,
            const Array *p_noise_y,
            glm::vec2    center,
            glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array               array = Array(shape);
  hmap::SlopeFunction f = hmap::SlopeFunction(angle, slope, center);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

Array smooth_cosine(glm::ivec2   shape,
                    const Array *p_noise_x,
                    const Array *p_noise_y,
                    glm::vec2    center,
                    glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array = Array(shape);

  auto lambda = [center](float x, float y, float)
  {
    float dx = x - center.x;
    float dy = y - center.y;
    float r = 2.f * M_PI * std::hypot(dx, dy);

    if (r < M_PI)
      return 0.5f + 0.5f * std::cos(r);
    else
      return 0.f;
  };

  fill_array_using_xy_function(array,
                               bbox,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               lambda);
  return array;
}

Array step(glm::ivec2   shape,
           float        angle,
           float        slope,
           const Array *p_ctrl_param,
           const Array *p_noise_x,
           const Array *p_noise_y,
           glm::vec2    center,
           glm::vec4    bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array              array = Array(shape);
  hmap::StepFunction f = hmap::StepFunction(angle, slope, center);

  fill_array_using_xy_function(array,
                               bbox,
                               p_ctrl_param,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               f.get_delegate());
  return array;
}

// --- Wrapper

Array get_primitive_base(const PrimitiveType &primitive_type,
                         const glm::ivec2    &shape,
                         const Array         *p_noise_x,
                         const Array         *p_noise_y,
                         glm::vec2            center,
                         glm::vec4            bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  switch (primitive_type)
  {
  case PrimitiveType::PRIM_BIQUAD_PULSE:
    return biquad_pulse(shape,
                        1.f,
                        nullptr,
                        p_noise_x,
                        p_noise_y,
                        center,
                        bbox);
    //
  case PrimitiveType::PRIM_BUMP:
    return bump(shape, 1.f, nullptr, p_noise_x, p_noise_y, center, bbox);
    //
  case PrimitiveType::PRIM_CONE:
    return cone(shape, 2.f, 1.f, false, center, p_noise_x, p_noise_y, bbox);
    //
  case PrimitiveType::PRIM_CONE_SMOOTH:
    return cone(shape, 2.f, 1.f, true, center, p_noise_x, p_noise_y, bbox);
  //
  case PrimitiveType::PRIM_CUBIC_PULSE:
    return cubic_pulse(shape, p_noise_x, p_noise_y, center, bbox);
    //
  case PrimitiveType::PRIM_SMOOTH_COSINE:
    return smooth_cosine(shape, p_noise_x, p_noise_y, center, bbox);
    //
  default: return Array(shape);
  }
}

} // namespace hmap
