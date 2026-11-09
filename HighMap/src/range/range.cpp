/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void chop(Array &array, float vmin)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [vmin](float x) { return x > vmin ? x : 0.f; };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void chop_max_smooth(Array &array, float vmax)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [vmax](float x)
  {
    if (x > vmax)
      x = 0.f;
    else if (x > 0.5f * vmax)
      x = vmax - x;
    return x;
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void clamp(Array &array, float vmin, float vmax)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [vmin, vmax](float x) { return std::clamp(x, vmin, vmax); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void clamp(Array &array, float vmax, ClampMode mode)
{
  if (!validate_non_empty(array)) return;

  switch (mode)
  {
  case ClampMode::POSITIVE_ONLY:
    // keep positives only, clamp max
    clamp(array, 0.f, vmax);
    return;

  case ClampMode::NEGATIVE_ONLY:
    // keep negatives only, reverse and clamp
    array *= -1.f;
    clamp(array, 0.f, vmax);
    return;

  case ClampMode::BOTH:
    // keep both, clamp symmetrically
    clamp(array, -vmax, vmax);
    return;

  case ClampMode::NONE:
  default: return;
  }
}

void clamp_max(Array &array, float vmax)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [vmax](float x) { return x < vmax ? x : vmax; };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void clamp_max(Array &array, const Array &vmax)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, vmax)) return;

  auto lambda = [](float x, float vmax) { return x < vmax ? x : vmax; };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 vmax.vector.begin(),
                 array.vector.begin(),
                 lambda);
}

void clamp_max_smooth(Array &array, float vmax, float k)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [k, vmax](float x)
  {
    float h = std::max(k - std::abs(x - vmax), 0.f) / k;
    return std::min(x, vmax) - h * h * h * k / 6.f;
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void clamp_max_smooth(Array &array, const Array &vmax, float k)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, vmax)) return;

  auto lambda = [k](float x, float vmax)
  {
    float h = std::max(k - std::abs(x - vmax), 0.f) / k;
    return std::min(x, vmax) - h * h * h * k / 6.f;
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 vmax.vector.begin(),
                 array.vector.begin(),
                 lambda);
}

void clamp_min(Array &array, float vmin)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [vmin](float x) { return x > vmin ? x : vmin; };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void clamp_min(Array &array, const Array &vmin)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, vmin)) return;

  auto lambda = [](float x, float vmin) { return x > vmin ? x : vmin; };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 vmin.vector.begin(),
                 array.vector.begin(),
                 lambda);
}

void clamp_min_smooth(Array &array, float vmin, float k)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [k, vmin](float x)
  {
    float h = std::max(k - std::abs(x - vmin), 0.f) / k;
    return std::max(x, vmin) + h * h * h * k / 6.f;
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void clamp_min_smooth(Array &array, const Array &vmin, float k)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, vmin)) return;

  auto lambda = [k](float x, float vmin)
  {
    float h = std::max(k - std::abs(x - vmin), 0.f) / k;
    return std::max(x, vmin) + h * h * h * k / 6.f;
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 vmin.vector.begin(),
                 array.vector.begin(),
                 lambda);
}

float clamp_min_smooth(float x, float vmin, float k)
{
  float h = std::max(k - std::abs(x - vmin), 0.f) / k;
  return std::max(x, vmin) + h * h * h * k / 6.f;
}

void clamp_oblique_plane(Array    &array,
                         float     vmax,
                         float     angle,
                         float     slope_value,
                         bool      use_max_operator,
                         float     k,
                         glm::vec2 center,
                         glm::vec4 bbox)
{
  if (!validate_non_empty(array)) return;

  // create plane
  Array plane = slope(array.shape,
                      angle,
                      slope_value,
                      /* p_ctrl_param  */ nullptr,
                      /* p_noise_x  */ nullptr,
                      /* p_noise_y  */ nullptr,
                      center,
                      bbox);
  plane += vmax;

  // clamp operation
  if (k == 0.f)
  {
    if (use_max_operator)
      clamp_max(array, plane);
    else
      clamp_min(array, plane);
  }
  else
  {
    if (use_max_operator)
      clamp_max_smooth(array, plane, k);
    else
      clamp_min_smooth(array, plane, k);
  }
}

void clamp_smooth(Array &array, float vmin, float vmax, float k)
{
  if (!validate_non_empty(array)) return;

  const float inv6k = k / 6.f;

  auto lambda = [k, vmin, vmax, inv6k](float x)
  {
    // smooth min
    float h = std::max(k - std::abs(x - vmin), 0.f) / k;
    float lo = std::max(x, vmin) + (h * h * h) * inv6k;

    // smooth max
    h = std::max(k - std::abs(lo - vmax), 0.f) / k;
    return std::min(lo, vmax) - (h * h * h) * inv6k;
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

Array maximum(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return std::max(a, b); });
  return array_out;
}

Array maximum(const Array &array1, const float value)
{
  if (!validate_non_empty(array1)) return Array();

  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array_out.vector.begin(),
                 [value](float a) { return std::max(a, value); });
  return array_out;
}

Array maximum_smooth(const Array &array1, const Array &array2, float k)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  if (k > 0.f)
  {
    Array array_out = Array(array1.shape);

    auto lambda = [k](float a, float b)
    {
      float h = std::max(k - std::abs(a - b), 0.f) / k;
      return std::max(a, b) + h * h * h * k / 6.f;
    };

    std::transform(array1.vector.begin(),
                   array1.vector.end(),
                   array2.vector.begin(),
                   array_out.vector.begin(),
                   lambda);
    return array_out;
  }
  else
    return maximum(array1, array2);
}

float maximum_smooth(const float a, const float b, float k)
{
  float h = std::max(k - std::abs(a - b), 0.f) / k;
  return std::max(a, b) + h * h * h * k / 6.f;
}

Array minimum(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array2.vector.begin(),
                 array_out.vector.begin(),
                 [](float a, float b) { return std::min(a, b); });
  return array_out;
}

Array minimum(const Array &array1, const float value)
{
  if (!validate_non_empty(array1)) return Array();

  Array array_out = Array(array1.shape);
  std::transform(array1.vector.begin(),
                 array1.vector.end(),
                 array_out.vector.begin(),
                 [value](float a) { return std::min(a, value); });
  return array_out;
}

Array minimum_smooth(const Array &array1, const Array &array2, float k)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  if (k > 0.f)
  {
    Array array_out = Array(array1.shape);

    auto lambda = [k](float a, float b)
    {
      float h = std::max(k - std::abs(a - b), 0.f) / k;
      return std::min(a, b) - h * h * h * k / 6.f;
    };

    std::transform(array1.vector.begin(),
                   array1.vector.end(),
                   array2.vector.begin(),
                   array_out.vector.begin(),
                   lambda);
    return array_out;
  }
  else
    return minimum(array1, array2);
}

float minimum_smooth(const float a, const float b, float k)
{
  float h = std::max(k - std::abs(a - b), 0.f) / k;
  return std::min(a, b) - h * h * h * k / 6.f;
}

void remap(Array &array, float vmin, float vmax, ICommunicator &comm)
{
  if (!validate_non_empty(array)) return;

  const glm::vec2 r = array.range(comm);
  const float     min = r.x;
  const float     max = r.y;

  if (min == max)
  {
    std::fill(array.vector.begin(), array.vector.end(), vmin);
    return;
  }

  const float scale = (vmax - vmin) / (max - min);
  const float offset = vmin - min * scale;

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 [scale, offset](float x) { return x * scale + offset; });
}

void remap(Array &array, float vmin, float vmax, float from_min, float from_max)
{
  if (!validate_non_empty(array)) return;

  if (from_min == from_max)
  {
    std::fill(array.vector.begin(), array.vector.end(), vmin);
    return;
  }

  const float scale = (vmax - vmin) / (from_max - from_min);
  const float offset = vmin - from_min * scale;

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 [scale, offset](float x) { return x * scale + offset; });
}

void rescale(Array &array, float scaling, float vref)
{
  if (!validate_non_empty(array)) return;

  if (vref == 0.f)
    // simply multiply the values by the scaling
    array *= scaling;
  else
  {
    const float offset = vref * (1.f - scaling);

    std::transform(array.vector.begin(),
                   array.vector.end(),
                   array.vector.begin(),
                   [scaling, offset](float x) { return x * scaling + offset; });
  }
}

void normalize_and_clip(Array         &array,
                        float          multiplier,
                        bool           smooth_transition,
                        bool           range_01,
                        ICommunicator &comm)
{
  if (!validate_non_empty(array)) return;

  float       mean = 0.f;
  const float sigma = array.std(&mean, comm);
  const float denom = std::max(multiplier * sigma, 1e-12f);

  if (range_01)
  {
    if (smooth_transition)
    {
      for (float &v : array.vector)
      {
        const float clamped = std::clamp((v - mean) / denom, -1.f, 1.f);
        const float t = 0.5f * (clamped + 1.f);

        v = smoothstep3(t);
      }
    }
    else
    {
      for (float &v : array.vector)
      {
        const float clamped = std::clamp((v - mean) / denom, -1.f, 1.f);

        v = 0.5f * (clamped + 1.f);
      }
    }
  }
  else
  {
    if (smooth_transition)
    {
      for (float &v : array.vector)
      {
        const float clamped = std::clamp((v - mean) / denom, -1.f, 1.f);
        const float t = 0.5f * (clamped + 1.f);

        v = 2.f * smoothstep3(t) - 1.f;
      }
    }
    else
    {
      for (float &v : array.vector)
        v = std::clamp((v - mean) / denom, -1.f, 1.f);
    }
  }
}

} // namespace hmap
