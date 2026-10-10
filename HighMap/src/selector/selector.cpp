/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/morphology.hpp"
#include "highmap/range.hpp"
#include "highmap/selector.hpp"

namespace hmap
{

Array scan_mask(const Array &array, float contrast, float brightness)
{
  if (!validate_non_empty(array)) return Array();

  brightness *= 0.5f;
  float low = std::clamp(contrast - brightness, 0.f, 1.f);
  float high = std::clamp(contrast + brightness, 0.f, 1.f);

  contrast = contrast * 2.f - 1.f;

  Array array_out = array + contrast;

  clamp(array_out);

  array_out -= low;
  array_out *= (high - low);

  array_out = smoothstep3(array_out);

  remap(array_out, 0.f, 1.f, 0.f, 1.f);

  return array_out;
}

Array select_angle(const Array &array, float angle, float sigma)
{
  if (!validate_non_empty(array)) return Array();
  if (sigma <= 0.f) return Array(array.shape, 0.f);

  Array           c = array;
  Array           aspect = gradient_angle(c, true);
  Array           gnorm = gradient_norm(c);
  constexpr float eps = 1e-9f;

  float target_rad = angle / 180.f * (float)M_PI;
  float sigma_rad = sigma / 180.f * (float)M_PI;

  Array res = Array(array.shape);
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      if (gnorm(i, j) <= eps)
      {
        res(i, j) = 0.f;
        continue;
      }

      float diff = std::atan2(std::sin(aspect(i, j) - target_rad),
                              std::cos(aspect(i, j) - target_rad));
      float d = std::abs(diff);
      res(i, j) = 1.f - threshold_smooth(d, 0.f, sigma_rad);
    }

  return res;
}

Array select_blob_log(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array c = array;
  smooth_cpulse(c, ir);
  c = -laplacian(c);
  extrapolate_borders(c, ir + 1, 0.1f);
  return c;
}

Array select_cavities(const Array &array, int ir, bool concave)
{
  if (!validate_non_empty(array)) return Array();

  Array array_smooth = array;
  smooth_cpulse(array_smooth, ir);
  Array c = -curvature_quadric(array_smooth, 0, CurvatureType::CT_MEAN);

  if (!concave) c *= -1.f;

  clamp_min(c, 0.f);
  return c;
}

Array select_elevation_slope(const Array &array,
                             float        gradient_scale,
                             float        vmax)
{
  if (!validate_non_empty(array)) return Array();

  Array da = gradient_norm(array) * gradient_scale * (float)array.shape.x;
  clamp_max(da, vmax);

  return pow((vmax - array) * (vmax - da), 0.5f);
}

Array select_elevation_slope(const Array &array, float gradient_scale)
{
  if (!validate_non_empty(array)) return Array();
  return select_elevation_slope(array, gradient_scale, array.max());
}

Array select_eq(const Array &array, float value)
{
  if (!validate_non_empty(array)) return Array();

  Array c = array;
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      c(i, j) = c(i, j) == value ? 1.f : 0.f;
  return c;
}

Array select_gt(const Array &array, float value)
{
  if (!validate_non_empty(array)) return Array();

  Array c = array;
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      c(i, j) = c(i, j) > value ? 1.f : 0.f;
  return c;
}

Array select_gradient_angle(const Array &array, float angle)
{
  if (!validate_non_empty(array)) return Array();

  Array c = gradient_angle(array);
  float alpha = angle / 180.f * M_PI;

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      c(i, j) = std::max(0.f, std::cos(alpha + c(i, j)));
  return c;
}

Array select_gradient_binary(const Array &array, float talus_center)
{
  if (!validate_non_empty(array)) return Array();

  Array c = gradient_norm(array);
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      c(i, j) = c(i, j) > talus_center ? 1.f : 0.f;
  return c;
}

Array select_gradient_exp(const Array &array,
                          float        talus_center,
                          float        talus_sigma)
{
  if (!validate_non_empty(array)) return Array();

  Array c = gradient_norm(array);
  c -= talus_center;
  c = exp(-c * c / (2.f * talus_sigma * talus_sigma));
  return c;
}

Array select_gradient_inv(const Array &array,
                          float        talus_center,
                          float        talus_sigma)
{
  if (!validate_non_empty(array)) return Array();

  Array c = gradient_norm(array);
  c -= talus_center;
  c = 1.f / (1.f + c * c / (talus_sigma * talus_sigma));
  return c;
}

Array select_interval(const Array &array, float value1, float value2)
{
  if (!validate_non_empty(array)) return Array();

  Array c = array;
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      if ((c(i, j) > value1) and (c(i, j) < value2))
        c(i, j) = 1.f;
      else
        c(i, j) = 0.f;
    }
  return c;
}

Array select_inward_outward_slope(const Array &array,
                                  glm::vec2    center,
                                  glm::vec4    bbox)
{
  if (!validate_non_empty(array)) return Array();

  Array c = Array(array.shape);

  glm::vec2 shift = {bbox.x, bbox.z};
  glm::vec2 scale = {bbox.y - bbox.x, bbox.w - bbox.z};

  int ic = (int)((center.x - shift.x) / scale.x * array.shape.x);
  int jc = (int)((center.y - shift.y) / scale.y * array.shape.y);

  for (int j = 0; j < array.shape.y - 1; j++)
    for (int i = 0; i < array.shape.x - 1; i++)
    {
      float hypot = (float)std::hypot(i - ic, j - jc);
      if (hypot > 0.f)
      {
        float u = (float)(i - ic) / hypot;
        float v = (float)(j - jc) / hypot;

        // elevation difference along the radial axis (if positive,
        // the slope is pointing to the center and is inward, otherwise the
        // slope is pointing outward)
        float dz = array.get_value_bilinear_at(i, j, u, v) - array(i, j);
        c(i, j) = dz;
      }
    }

  extrapolate_borders(c);

  return c;
}

Array select_lt(const Array &array, float value)
{
  if (!validate_non_empty(array)) return Array();

  Array c = array;
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      c(i, j) = c(i, j) < value ? 1.f : 0.f;
  return c;
}

Array select_midrange(const Array &array, float gain, float vmin, float vmax)
{
  if (!validate_non_empty(array)) return Array();

  Array c = array;
  hmap::remap(c, -1.f, 1.f, vmin, vmax);

  float norm_coeff = 1.f / std::exp(-1.f);

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      float v = c(i, j) * c(i, j);
      if (v != 1.f)
        c(i, j) = std::pow(norm_coeff * std::exp(-1.f / (1.f - v)), 1.f / gain);
      else
        c(i, j) = 0.f;
    }

  return c;
}

Array select_midrange(const Array &array, float gain)
{
  if (!validate_non_empty(array)) return Array();
  return select_midrange(array, gain, array.min(), array.max());
}

void select_multiband3(const Array &array,
                       Array       &band_low,
                       Array       &band_mid,
                       Array       &band_high,
                       float        ratio1,
                       float        ratio2,
                       float        overlap,
                       float        vmin,
                       float        vmax)
{
  if (!validate_non_empty(array)) return;

  band_low = Array(array.shape);
  band_mid = Array(array.shape);
  band_high = Array(array.shape);

  float v1 = vmin + ratio1 * (vmax - vmin);
  float v2 = vmin + ratio2 * (vmax - vmin);

  auto lambda = [](float r, float r0, float r1, float r2, float overlap)
  {
    float w0 = overlap * (r1 - r0);
    float w2 = overlap * (r2 - r1);

    if (r < r0 - w0) return 0.f;
    if (r > r2 + w2) return 0.f;

    if (r > r0 + w0 && r < r2 - w2) return 1.f;

    float rn = 0.f;
    if (r < r0 + w0) rn = (r - r0 + w0) / 2.f / w0;
    if (r > r2 - w2) rn = 1.f - (r - r2 + w2) / 2.f / w2;
    return rn * rn * (3.f - 2.f * rn);
  };

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      band_low(i,
               j) = lambda(array(i, j), vmin, 0.5f * (vmin + v1), v1, overlap);
      band_mid(i, j) = lambda(array(i, j), v1, 0.5f * (v1 + v2), v2, overlap);
      band_high(i,
                j) = lambda(array(i, j), v2, 0.5f * (v2 + vmax), vmax, overlap);
    }
}

void select_multiband3(const Array &array,
                       Array       &band_low,
                       Array       &band_mid,
                       Array       &band_high,
                       float        ratio1,
                       float        ratio2,
                       float        overlap)
{
  if (!validate_non_empty(array)) return;

  float vmin = array.min();
  float vmax = array.max();

  select_multiband3(array,
                    band_low,
                    band_mid,
                    band_high,
                    ratio1,
                    ratio2,
                    overlap,
                    vmin,
                    vmax);
}

Array select_pulse(const Array &array, float value, float sigma)
{
  if (!validate_non_empty(array)) return Array();

  Array c = Array(array.shape);
  float a = 1.f / sigma;
  float b = -value / sigma;

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      float r = std::abs(a * array(i, j) + b);
      if (r < 1.f) c(i, j) = 1.f - r * r * (3.f - 2.f * r);
    }
  return c;
}

Array select_range(const Array &array, float vmin, float vmax, float width)
{
  if (!validate_non_empty(array)) return Array();

  float low = std::min(vmin, vmax);
  float high = std::max(vmin, vmax);
  float w = std::max(0.0f, width);

  Array out(array.shape);

  if (w <= 1e-7f)
  {
    for (int j = 0; j < array.shape.y; ++j)
      for (int i = 0; i < array.shape.x; ++i)
      {
        float v = array(i, j);
        out(i, j) = (v >= low && v <= high) ? 1.0f : 0.0f;
      }
    return out;
  }

  for (int j = 0; j < array.shape.y; ++j)
    for (int i = 0; i < array.shape.x; ++i)
    {
      float v = array(i, j);
      float left = threshold_smooth(v, low - w, low);
      float right = 1.0f - threshold_smooth(v, high, high + w);
      out(i, j) = left * right;
    }

  return out;
}

Array select_rivers(const Array &array, float talus_ref, float clipping_ratio)
{
  if (!validate_non_empty(array)) return Array();

  // see erosion/hydraulic_stream
  Array facc = flow_accumulation_dinf(array, talus_ref);
  float vmax = clipping_ratio * std::pow(facc.sum() / (float)facc.size(), 0.5f);

  if (vmax == 0.f) // failsafe
    return Array(array.shape);

  clamp(facc, 0.f, vmax);
  return facc / vmax;
}

Array select_transitions(const Array &array1,
                         const Array &array2,
                         const Array &array_blend)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2) ||
      !validate_same_shape(array1, array_blend))
    return Array();

  // set the whole mask to 1 and look for "non-transitioning" regions
  Array mask = Array(array1.shape, 1.f);

  for (int j = 0; j < array1.shape.y - 1; j++)
    for (int i = 0; i < array1.shape.x - 1; i++)
    {
      if ((array_blend(i, j) == array1(i, j)) &
          (array_blend(i + 1, j) == array1(i + 1, j)) &
          (array_blend(i, j) == array1(i, j)) &
          (array_blend(i, j + 1) == array1(i, j + 1)))
        mask(i, j) = 0.f;

      else if ((array_blend(i, j) == array2(i, j)) &
               (array_blend(i + 1, j) == array2(i + 1, j)) &
               (array_blend(i, j) == array2(i, j)) &
               (array_blend(i, j + 1) == array2(i, j + 1)))
        mask(i, j) = 0.f;
    }

  // boundaries
  for (int j = 0; j < array1.shape.y - 1; j++)
  {
    int i = array1.shape.x - 1;

    if ((array_blend(i, j) == array1(i, j)) &
        (array_blend(i - 1, j) == array1(i - 1, j)) &
        (array_blend(i, j) == array1(i, j)) &
        (array_blend(i, j + 1) == array1(i, j + 1)))
      mask(i, j) = 0.f;

    else if ((array_blend(i, j) == array2(i, j)) &
             (array_blend(i - 1, j) == array2(i - 1, j)) &
             (array_blend(i, j) == array2(i, j)) &
             (array_blend(i, j + 1) == array2(i, j + 1)))
      mask(i, j) = 0.f;
  }

  for (int i = 0; i < array1.shape.x - 1; i++)
  {
    int j = array1.shape.y - 1;

    if ((array_blend(i, j) == array1(i, j)) &
        (array_blend(i + 1, j) == array1(i + 1, j)) &
        (array_blend(i, j) == array1(i, j)) &
        (array_blend(i, j - 1) == array1(i, j - 1)))
      mask(i, j) = 0.f;

    else if ((array_blend(i, j) == array2(i, j)) &
             (array_blend(i + 1, j) == array2(i + 1, j)) &
             (array_blend(i, j) == array2(i, j)) &
             (array_blend(i, j - 1) == array2(i, j - 1)))
      mask(i, j) = 0.f;
  }

  return mask;
}

Array select_valley(const Array &z, int ir, bool ridge_select)
{
  if (!validate_non_empty(z)) return Array();

  if (ridge_select)
    return morphological_top_hat(z, ir);
  else
    return morphological_black_hat(z, ir);
}

Array smooth_mask_preserve_frontier(const Array &mask, float feather_width)
{
  if (!validate_non_empty(mask)) return Array();

  if (feather_width <= 0.f) return mask;

  // distance to the 0/1 frontier inside the mask
  Array dist = distance_transform(is_zero(mask));

  Array smoothed = mask;

  for (int j = 0; j < mask.shape.y; j++)
    for (int i = 0; i < mask.shape.x; i++)
    {
      if (mask(i, j) > 0.f)
      {
        float t = std::clamp(dist(i, j) / feather_width, 0.f, 1.f);
        smoothed(i, j) = smoothstep3(t);
      }
    }

  return smoothed;
}

} // namespace hmap
