/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"

#include <format>

namespace hmap
{

Array abs(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::abs(v); });
  return array_out;
}

Array abs_smooth(const Array &array, float k)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  float k2 = k * k;
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [&k2](float v) { return std::sqrt(v * v + k2); });
  return array_out;
}

Array abs_smooth(const Array &array, float k, float vshift)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  float k2 = k * k;
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [&k2, &vshift](float v)
                 {
                   float vbis = v - vshift;
                   return vshift + std::sqrt(vbis * vbis + k2);
                 });
  return array_out;
}

Array abs_smooth(const Array &array, float k, const Array &vshift)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_same_shape(array, vshift)) return Array();

  Array array_out = Array(array.shape);
  float k2 = k * k;
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 vshift.vector.begin(),
                 array_out.vector.begin(),
                 [&k2](float v, float s)
                 {
                   float vbis = v - s;
                   return s + std::sqrt(vbis * vbis + k2);
                 });
  return array_out;
}

Array almost_unit_identity(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return almost_unit_identity(v); });
  return array_out;
}

Array almost_unit_identity_c2(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return almost_unit_identity_c2(v); });
  return array_out;
}

Array atan(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::atan(v); });
  return array_out;
}

Array atan2(const Array &y, const Array &array)
{
  if (!validate_non_empty(y)) return Array();
  if (!validate_same_shape(y, array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(y.vector.begin(),
                 y.vector.end(),
                 array.vector.begin(),
                 array_out.vector.begin(),
                 [](float y_, float x_) { return std::atan2(y_, x_); });
  return array_out;
}

Array build_density_linear(const std::vector<const Array *> &layers,
                           const std::vector<float>         &weights)
{
  if (layers.empty()) return Array();
  if (layers[0] == nullptr || !validate_non_empty(*layers[0])) return Array();
  if (!weights.empty() && weights.size() != layers.size()) return Array();

  for (size_t k = 1; k < layers.size(); ++k)
  {
    if (layers[k] == nullptr || !validate_same_shape(*layers[0], *layers[k]))
      return Array();
  }

  const glm::ivec2 shape = layers[0]->shape;
  Array            array_out(shape, 1.f);
  const size_t     n = layers[0]->vector.size();

  for (size_t k = 0; k < layers.size(); ++k)
  {
    float w = weights.empty() ? 1.f : weights[k];
    if (w == 0.f) continue;

    const auto &vec = layers[k]->vector;
    if (w == 1.f)
    {
      for (size_t i = 0; i < n; ++i)
        array_out.vector[i] *= vec[i];
    }
    else
    {
      for (size_t i = 0; i < n; ++i)
        array_out.vector[i] *= std::lerp(1.f, vec[i], w);
    }
  }

  return array_out;
}

Array build_density_log(const std::vector<const Array *> &layers,
                        const std::vector<float>         &weights,
                        float                             eps)
{
  if (layers.empty()) return Array();
  if (layers[0] == nullptr || !validate_non_empty(*layers[0])) return Array();
  if (!weights.empty() && weights.size() != layers.size()) return Array();

  for (size_t k = 1; k < layers.size(); ++k)
  {
    if (layers[k] == nullptr || !validate_same_shape(*layers[0], *layers[k]))
      return Array();
  }

  const glm::ivec2   shape = layers[0]->shape;
  const size_t       n = layers[0]->vector.size();
  std::vector<float> sum_log(n, 0.f);

  for (size_t k = 0; k < layers.size(); ++k)
  {
    float w = weights.empty() ? 1.f : weights[k];
    if (w == 0.f) continue;

    const auto &vec = layers[k]->vector;
    for (size_t i = 0; i < n; ++i)
      sum_log[i] += w * log_safe(vec[i], eps);
  }

  Array array_out(shape);
  for (size_t i = 0; i < n; ++i)
    array_out.vector[i] = std::exp(sum_log[i]);

  return array_out;
}

Array cos(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::cos(v); });
  return array_out;
}

Array exp(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::exp(v); });
  return array_out;
}

Array gaussian_decay(const Array &array, float sigma)
{
  if (!validate_non_empty(array)) return Array();

  Array       array_out = Array(array.shape);
  const float coeff = 0.5f / (sigma * sigma);

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [coeff](float v) { return std::exp(-coeff * v * v); });
  return array_out;
}

Array hypot(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1)) return Array();
  if (!validate_same_shape(array1, array2)) return Array();

  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return std::sqrt(a * a + b * b); });
  return array_out;
}

Array is_equal(const Array &array, float value)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [value](float v) { return v == value ? 1.f : 0.f; });
  return array_out;
}

Array is_non_zero(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return v != 0.f ? 1.f : 0.f; });
  return array_out;
}

Array is_zero(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return v == 0.f ? 1.f : 0.f; });
  return array_out;
}

Array lerp(const Array &array1, const Array &array2, const Array &t)
{
  if (!validate_non_empty(array1)) return Array();
  if (!validate_same_shape(array1, array2)) return Array();
  if (!validate_same_shape(array1, t)) return Array();

  Array        array_out(array1.shape);
  const size_t n = array1.vector.size();

  for (size_t i = 0; i < n; ++i)
  {
    float ti = t.vector[i];
    array_out.vector[i] = array1.vector[i] * (1.f - ti) + array2.vector[i] * ti;
  }

  return array_out;
}

Array lerp(const Array &array1, const Array &array2, float t)
{
  if (!validate_non_empty(array1)) return Array();
  if (!validate_same_shape(array1, array2)) return Array();

  Array        array_out(array1.shape);
  const size_t n = array1.vector.size();

  for (size_t i = 0; i < n; ++i)
    array_out.vector[i] = array1.vector[i] * (1.f - t) + array2.vector[i] * t;

  return array_out;
}

Array log10(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::log10(v); });
  return array_out;
}

Array log_safe(const Array &array, float eps)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [eps](float v) { return log_safe(v, eps); });
  return array_out;
}

Array pow(const Array &array, float exp)
{
  if (!validate_non_empty(array)) return Array();

  if (exp == 2.f) return array * array;
  if (exp == 0.5f) return sqrt(array);
  if (exp == 3.f)
  {
    Array sq = array * array;
    return sq * array;
  }

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [exp](float v) { return std::pow(v, exp); });
  return array_out;
}

Array r_min(const Array &array1, const Array &array2, float alpha)
{
  if (!validate_non_empty(array1)) return Array();
  if (!validate_same_shape(array1, array2)) return Array();

  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 [alpha](float a, float b) { return r_min(a, b, alpha); });
  return array_out;
}

Array r_max(const Array &array1, const Array &array2, float alpha)
{
  if (!validate_non_empty(array1)) return Array();
  if (!validate_same_shape(array1, array2)) return Array();

  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 [alpha](float a, float b) { return r_max(a, b, alpha); });
  return array_out;
}

Array sigmoid(const Array &array, float width, float vmin, float vmax, float x0)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [width, vmin, vmax, x0](float x)
                 { return sigmoid(x, width, vmin, vmax, x0); });
  return array_out;
}

Array sin(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::sin(v); });
  return array_out;
}

Array smoothstep3(const Array &array, float vmin, float vmax)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);

  const float inv_range = 1.f / (vmax - vmin);
  const float range = vmax - vmin;

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [vmin, vmax, range, inv_range](float v)
                 {
                   if (v <= vmin) return vmin;
                   if (v >= vmax) return vmax;
                   float vn = (v - vmin) * inv_range;
                   vn = vn * vn * (3.f - 2.f * vn);
                   return vmin + range * vn;
                 });

  return array_out;
}

Array smoothstep3_lower(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return smoothstep3_lower(v); });
  return array_out;
}

Array smoothstep3_upper(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return smoothstep3_upper(v); });
  return array_out;
}

Array smoothstep5(const Array &array, float vmin, float vmax)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [vmin, vmax](float v)
                 {
                   if (v < vmin)
                     return vmin;
                   else if (v > vmax)
                     return vmax;
                   else
                   {
                     float vn = (v - vmin) / (vmax - vmin);
                     vn = vn * vn * vn * (vn * (vn * 6.f - 15.f) + 10.f);
                     return vmin + (vmax - vmin) * vn;
                   }
                 });
  return array_out;
}

Array smoothstep5(const Array &array, const Array &vmin, const Array &vmax)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_same_shape(array, vmin)) return Array();
  if (!validate_same_shape(array, vmax)) return Array();

  Array        array_out(array.shape);
  const size_t n = array.vector.size();

  for (size_t i = 0; i < n; ++i)
  {
    const float v = array.vector[i];
    const float lo = vmin.vector[i];
    const float hi = vmax.vector[i];

    if (v <= lo)
    {
      array_out.vector[i] = lo;
      continue;
    }
    else if (v >= hi)
    {
      array_out.vector[i] = hi;
      continue;
    }

    float vn = (v - lo) / (hi - lo);
    vn = vn * vn * vn * (vn * (vn * 6.f - 15.f) + 10.f);
    array_out.vector[i] = lo + (hi - lo) * vn;
  }
  return array_out;
}

Array smoothstep5_lower(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return smoothstep5_lower(v); });
  return array_out;
}

Array smoothstep5_upper(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return smoothstep5_upper(v); });
  return array_out;
}

Array smoothstep7(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return smoothstep7(v); });
  return array_out;
}

Array sqrt(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return std::sqrt(v); });
  return array_out;
}

Array sqrt_safe(const Array &array)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [](float v) { return v > 0.f ? std::sqrt(v) : 0.f; });
  return array_out;
}

Array threshold(const Array &array, float x0, float x1)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [x0, x1](float v) { return threshold(v, x0, x1); });
  return array_out;
}

Array threshold_smooth(const Array &array, float x0, float x1)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);
  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array_out.vector.begin(),
                 [x0, x1](float v) { return threshold_smooth(v, x0, x1); });
  return array_out;
}

} // namespace hmap
