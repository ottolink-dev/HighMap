/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "highmap/math/core.hpp"

namespace hmap
{

float abs_smooth(float a, float k)
{
  float k2 = k * k;
  return std::sqrt(a * a + k2);
}

float almost_unit_identity(float x)
{
  return (2.f - x) * x * x;
}

float almost_unit_identity_c2(float x)
{
  // second-order derivative equals 0 at x = 1 also to avoid
  // discontinuities in some cases
  return x * x * (x * x - 3.f * x + 3.f);
}

int ceil_div(int a, int b)
{
  return (a + b - 1) / b;
}

float fast_exp(float x)
{
  union
  {
    float   f;
    int32_t i;
  } u;
  u.i = (int32_t)(12102203.f * x) + 127 * (1 << 23);
  return u.f;
}

float fast_log(float x)
{
  union
  {
    float   f;
    int32_t i;
  } u = {x};
  return (u.i - 1064866805) * 8.262958405e-8f;
}

float gain(float x, float factor)
{
  return x < 0.5 ? 0.5f * std::pow(2.f * x, factor)
                 : 1.f - 0.5f * std::pow(2.f * (1.f - x), factor);
}

float lerp(float a, float b, float t)
{
  return std::lerp(a, b, t);
}

float power_curve(float x, float a, float b)
{
  // https://iquilezles.org/articles/functions/
  float k = std::pow(a + b, a + b) / (std::pow(a, a) * std::pow(b, b));
  return k * std::pow(x, a) * std::pow(1.f - x, b);
}

float r_min(float f1, float f2, float alpha)
{
  float disc = std::max(0.f, f1 * f1 + f2 * f2 - 2.f * alpha * f1 * f2);
  return (f1 + f2 - std::sqrt(disc)) / (1.f + alpha);
}

float r_max(float f1, float f2, float alpha)
{
  float disc = std::max(0.f, f1 * f1 + f2 * f2 - 2.f * alpha * f1 * f2);
  return (f1 + f2 + std::sqrt(disc)) / (1.f + alpha);
}

float sigmoid(float x, float width, float vmin, float vmax, float x0)
{
  float v = 1.f / (1.f + std::exp(-(x - x0) / width));
  v = (vmax - vmin) * v + vmin;
  return v;
}

float smoothstep3(float x)
{
  return x * x * (3.f - 2.f * x);
}

float smoothstep3_lower(float x)
{
  return x * (2.f * x - x * x);
}

float smoothstep3_upper(float x)
{
  return x * (1.f + x - x * x);
}

float smoothstep5(float x)
{
  return x * x * x * (x * (x * 6.f - 15.f) + 10.f);
}

float smoothstep5_lower(float x)
{
  return x * x * x * (6.f - 8.f * x + 3.f * x * x);
}

float smoothstep5_upper(float x)
{
  return x * (1.f + x * x * (4.f - 7.f * x + 3.f * x * x));
}

float smoothstep7(float x)
{
  float x2 = x * x;
  float x3 = x2 * x;
  float x4 = x3 * x;
  float x5 = x4 * x;
  float x6 = x5 * x;
  float x7 = x6 * x;
  return -20.f * x7 + 70.f * x6 - 84.f * x5 + 35.f * x4;
}

float threshold(float x, float x0, float x1)
{
  if (x < x0) return 0.f;
  if (x < x1) return (x - x0) / (x1 - x0);
  return 1.f;
}

float threshold_smooth(float x, float x0, float x1)
{
  if (x < x0) return 0.f;
  if (x < x1) return smoothstep3((x - x0) / (x1 - x0));
  return 1.f;
}

float trapeze(float x, float x0, float x1, float width)
{
  if (width <= 0.0f) return (x >= x0 && x <= x1) ? 1.0f : 0.0f;

  // rising edge
  float left = (x - (x0 - width)) / width;

  // falling edge
  float right = ((x1 + width) - x) / width;

  return std::clamp(std::min(left, right), 0.0f, 1.0f);
}

float trapeze_smooth(float x, float x0, float x1, float width)
{
  return smoothstep3(trapeze(x, x0, x1, width));
}

float triangle(float x, float vmin, float vmax)
{
  if (x <= vmin || x >= vmax) return 0.f;

  float mid = 0.5f * (vmin + vmax);
  float half_width = 0.5f * (vmax - vmin);

  return 1.f - std::abs((x - mid) / half_width);
}

} // namespace hmap
