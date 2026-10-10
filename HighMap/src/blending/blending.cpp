/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/range.hpp"

#include <format>

namespace hmap
{

Array blend_exclusion(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array array_out = Array(array1.shape);
  array_out = 0.5f - 2.f * (-0.5f + array1) * (-0.5f + array2);
  return array_out;
}

Array blend_gradients(const Array &array1, const Array &array2, int ir)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array dn1 = gradient_norm(array1);
  Array dn2 = gradient_norm(array2);

  smooth_cpulse(dn1, ir);
  smooth_cpulse(dn2, ir);

  Array t = maximum_smooth(dn1, dn2, 0.1f);
  remap(t);

  return lerp(array1, array2, t);
}

Array blend_negate(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array array_out = Array(array1.shape);

  auto lambda = [](float a, float b) { return a < b ? a : 2.f * b - a; };

  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 lambda);

  return array_out;
}

Array blend_overlay(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array array_out = Array(array1.shape);
  auto  lambda = [](float a, float b)
  { return a < 0.5 ? 2.f * a * b : 1.f - 2.f * (1.f - a) * (1.f - b); };

  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 lambda);

  return array_out;
}

Array blend_soft(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array array_out = Array(array1.shape);
  array_out = (1.f - array1) * array1 * array2 +
              array1 * (1.f - (1.f - array1) * (1.f - array2));
  return array_out;
}

Array blend_power_law(const std::vector<const Array *> &arrays, float alpha)
{
  if (arrays.empty()) return Array();
  if (arrays[0] == nullptr || !validate_non_empty(*arrays[0])) return Array();

  for (size_t k = 1; k < arrays.size(); ++k)
  {
    if (arrays[k] == nullptr || !validate_same_shape(*arrays[0], *arrays[k]))
      return Array();
  }

  if (arrays.size() == 1) return *arrays[0];

  glm::ivec2 shape = arrays[0]->shape;
  Array      num = Array(shape, 0.f);
  Array      den = Array(shape, 0.f);

  if (alpha == 0.f)
  {
    for (const auto *arr : arrays)
      num += *arr;
    return num / static_cast<float>(arrays.size());
  }

  for (const auto *arr : arrays)
  {
    Array p = pow(*arr, alpha);
    den += p;
    num += p * (*arr);
  }

  return num / den;
}

Array blend_power_law(const Array &array1, const Array &array2, float alpha)
{
  return blend_power_law({&array1, &array2}, alpha);
}

Array mixer(const Array                      &t,
            const std::vector<const Array *> &arrays,
            float                             gain_factor)
{
  if (!validate_non_empty(t)) return Array();
  if (arrays.empty()) return Array(t.shape);

  for (const auto *arr : arrays)
  {
    if (arr == nullptr || !validate_same_shape(t, *arr)) return Array(t.shape);
  }

  if (arrays.size() == 1) return *arrays[0];

  Array               array_out = Array(t.shape);
  const std::uint32_t n = arrays.size();

  for (std::uint32_t k = 0; k < n; k++)
  {
    float r0 = (float)k / (float)(n - 1);

    if (gain_factor == 1.f)
    {
      for (int j = 0; j < t.shape.y; j++)
        for (int i = 0; i < t.shape.x; i++)
        {
          float ta = 1.f - std::fabs(t(i, j) - r0) * (float)(n - 1);
          if (ta >= 0.f)
          {
            float ts = ta * ta * (3.f - 2.f * ta);
            array_out(i, j) += ts * (*arrays[k])(i, j);
          }
        }
    }
    else
    {
      for (int j = 0; j < t.shape.y; j++)
        for (int i = 0; i < t.shape.x; i++)
        {
          float ta = 1.f - std::fabs(t(i, j) - r0) * (float)(n - 1);
          if (ta >= 0.f)
          {
            float ts = ta * ta * (3.f - 2.f * ta);
            ts = gain(ts, gain_factor);
            array_out(i, j) += ts * (*arrays[k])(i, j);
          }
        }
    }
  }
  return array_out;
}

Array transfer(const Array &source,
               const Array &target,
               int          ir,
               float        amplitude,
               bool         target_prefiltering)
{
  if (!validate_non_empty(source) || !validate_same_shape(source, target))
    return Array();

  // high-pass spatial filter
  Array w = -source;
  smooth_cpulse(w, ir);
  w += source;

  if (target_prefiltering)
  {
    Array target_f = target;
    smooth_cpulse(target_f, ir);
    w = target_f + amplitude * w;
  }
  else
  {
    w = target + amplitude * w;
  }

  return w;
}

} // namespace hmap
