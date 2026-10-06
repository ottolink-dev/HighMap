/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/interpolate/interpolate_array.hpp"
#include "highmap/logger.hpp"
#include "highmap/operator.hpp"
#include "highmap/random.hpp"

namespace hmap
{

// NB - Interpolating from an empty source has no defined result, and the index
// arithmetic in these routines is undefined for a zero dimension:
// std::clamp(v, 0, shape.x - 1) becomes std::clamp(v, 0, -1), whose contract is
// violated when lo > hi and which in practice returns -1, indexing before the
// buffer. 1 / shape.x is also inf. Leave the target as-is — it is zero-filled
// by its constructor, which is the only sensible result here.

void interpolate_array_bicubic(const Array &source,
                               Array       &target,
                               bool         endpoint,
                               bool         pixel_centered)
{
  glm::vec4 bbox_source(0.f, 1.f, 0.f, 1.f);
  glm::vec4 bbox_target(0.f, 1.f, 0.f, 1.f);

  interpolate_array_bicubic(source,
                            target,
                            bbox_source,
                            bbox_target,
                            endpoint,
                            pixel_centered);
}

void interpolate_array_bicubic(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target,
                               bool             endpoint,
                               bool             pixel_centered)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  float dx_s = 1.f / static_cast<float>(source.shape.x);
  float dy_s = 1.f / static_cast<float>(source.shape.y);

  float dx_t = 1.f / static_cast<float>(target.shape.x);
  float dy_t = 1.f / static_cast<float>(target.shape.y);

  float shift_x_t = pixel_centered ? 0.5f * dx_t : 0.f;
  float shift_y_t = pixel_centered ? 0.5f * dy_t : 0.f;
  float shift_x_s = pixel_centered ? 0.5f : 0.f;
  float shift_y_s = pixel_centered ? 0.5f : 0.f;

  // grid points (pixel-centered)
  std::vector<float> x = linspace(bbox_target.x + shift_x_t,
                                  bbox_target.y,
                                  target.shape.x,
                                  endpoint);

  std::vector<float> y = linspace(bbox_target.z + shift_y_t,
                                  bbox_target.w,
                                  target.shape.y,
                                  endpoint);

  // rescale to a unit square based on bbox_source
  for (auto &x_ : x)
    x_ = (x_ - bbox_source.x) / (bbox_source.y - bbox_source.x);

  for (auto &y_ : y)
    y_ = (y_ - bbox_source.z) / (bbox_source.w - bbox_source.z);

  for (int j = 0; j < target.shape.y; ++j)
    for (int i = 0; i < target.shape.x; ++i)
    {
      // reference source index
      float xc = x[i] / dx_s - shift_x_s;
      float yc = y[j] / dy_s - shift_y_s;

      // Nb - no clamping because it will be lead to issues with u & v
      // calculations and clamping is actually done later on while
      // interpolating
      int is0 = static_cast<int>(xc);
      int js0 = static_cast<int>(yc);

      // interp
      float u = xc - is0;
      float v = yc - js0;

      // interpolate
      float arr[4][4];

      // get the 4x4 surrounding grid points
      for (int n = -1; n <= 2; ++n)
        for (int m = -1; m <= 2; ++m)
        {
          int ip = std::clamp(is0 + m, 0, source.shape.x - 1);
          int jp = std::clamp(js0 + n, 0, source.shape.y - 1);
          arr[m + 1][n + 1] = source(ip, jp);
        }

      // interpolate in the x direction
      float col_results[4];
      for (int k = 0; k < 4; ++k)
        col_results[k] = cubic_interpolate(arr[k], v);

      // interpolate in the y direction
      target(i, j) = cubic_interpolate(col_results, u);
    }
}

void interpolate_array_bilinear(const Array &source,
                                Array       &target,
                                bool         endpoint,
                                bool         pixel_centered)
{
  glm::vec4 bbox_source(0.f, 1.f, 0.f, 1.f);
  glm::vec4 bbox_target(0.f, 1.f, 0.f, 1.f);

  interpolate_array_bilinear(source,
                             target,
                             bbox_source,
                             bbox_target,
                             endpoint,
                             pixel_centered);
}

void interpolate_array_bilinear(const Array     &source,
                                Array           &target,
                                const glm::vec4 &bbox_source,
                                const glm::vec4 &bbox_target,
                                bool             endpoint,
                                bool             pixel_centered)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  float dx_s = 1.f / static_cast<float>(source.shape.x);
  float dy_s = 1.f / static_cast<float>(source.shape.y);

  float dx_t = 1.f / static_cast<float>(target.shape.x);
  float dy_t = 1.f / static_cast<float>(target.shape.y);

  float shift_x_t = pixel_centered ? 0.5f * dx_t : 0.f;
  float shift_y_t = pixel_centered ? 0.5f * dy_t : 0.f;
  float shift_x_s = pixel_centered ? 0.5f : 0.f;
  float shift_y_s = pixel_centered ? 0.5f : 0.f;

  // grid points (pixel-centered)
  std::vector<float> x = linspace(bbox_target.x + shift_x_t,
                                  bbox_target.y,
                                  target.shape.x,
                                  endpoint);

  std::vector<float> y = linspace(bbox_target.z + shift_y_t,
                                  bbox_target.w,
                                  target.shape.y,
                                  endpoint);

  // rescale to a unit square based on bbox_source
  for (auto &x_ : x)
    x_ = (x_ - bbox_source.x) / (bbox_source.y - bbox_source.x);

  for (auto &y_ : y)
    y_ = (y_ - bbox_source.z) / (bbox_source.w - bbox_source.z);

  for (int j = 0; j < target.shape.y; ++j)
    for (int i = 0; i < target.shape.x; ++i)
    {
      // reference source index
      float xc = x[i] / dx_s - shift_x_s;
      float yc = y[j] / dy_s - shift_y_s;

      int is0 = static_cast<int>(xc);
      int js0 = static_cast<int>(yc);

      is0 = std::clamp(is0, 0, source.shape.x - 1);
      js0 = std::clamp(js0, 0, source.shape.y - 1);

      // interp
      float u = xc - is0;
      float v = yc - js0;

      int is1 = std::min(is0 + 1, source.shape.x - 1);
      int js1 = std::min(js0 + 1, source.shape.y - 1);

      target(i, j) = bilinear_interp(source(is0, js0),
                                     source(is1, js0),
                                     source(is0, js1),
                                     source(is1, js1),
                                     u,
                                     v);
    }
}

void interpolate_array_nearest(const Array &source,
                               Array       &target,
                               bool         endpoint)
{
  glm::vec4 bbox_source(0.f, 1.f, 0.f, 1.f);
  glm::vec4 bbox_target(0.f, 1.f, 0.f, 1.f);

  interpolate_array_nearest(source, target, bbox_source, bbox_target, endpoint);
}

void interpolate_array_nearest(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target,
                               bool             endpoint)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  std::vector<float> x = linspace(bbox_target.x,
                                  bbox_target.y,
                                  target.shape.x,
                                  endpoint);

  std::vector<float> y = linspace(bbox_target.z,
                                  bbox_target.w,
                                  target.shape.y,
                                  endpoint);

  for (int j = 0; j < target.shape.y; ++j)
    for (int i = 0; i < target.shape.x; ++i)
    {
      int is = static_cast<int>(
          std::round((x[i] - bbox_source.x) / (bbox_source.y - bbox_source.x) *
                     source.shape.x));
      int js = static_cast<int>(
          std::round((y[j] - bbox_source.z) / (bbox_source.w - bbox_source.z) *
                     source.shape.y));

      is = std::clamp(is, 0, source.shape.x - 1);
      js = std::clamp(js, 0, source.shape.y - 1);

      target(i, j) = source(is, js);
    }
}

Array resample_bicubic_warp(const Array  &src,
                            glm::ivec2    dst_shape,
                            float         warp,
                            std::uint32_t seed)
{
  if (!validate_non_empty(src) || !validate_shape(dst_shape)) return Array();

  if (warp <= 0.f) return src.resample_to_shape_bicubic(dst_shape);

  int sw = src.shape.x;
  int sh = src.shape.y;
  int dw = dst_shape.x;
  int dh = dst_shape.y;

  Array dst(dst_shape);

  const float rx = static_cast<float>(sw) / static_cast<float>(dw);
  const float ry = static_cast<float>(sh) / static_cast<float>(dh);

#pragma omp parallel for collapse(2) schedule(static)
  for (int y = 0; y < dh; ++y)
  {
    for (int x = 0; x < dw; ++x)
    {
      float fx = (static_cast<float>(x) + 0.5f) * rx - 0.5f;
      float fy = (static_cast<float>(y) + 0.5f) * ry - 0.5f;

      if (warp > 0.f)
      {
        float u = fx * 0.5f;
        float v = fy * 0.5f;
        fx += warp * vnoise(u + 13.7f, v + 3.1f, seed);
        fy += warp * vnoise(u + 71.3f, v + 47.9f, seed + 101u);
      }

      int   xi = static_cast<int>(std::floor(fx));
      int   yi = static_cast<int>(std::floor(fy));
      float wx[4], wy[4];
      bspline_weights(fx - static_cast<float>(xi), wx);
      bspline_weights(fy - static_cast<float>(yi), wy);

      float val = 0.f;
      for (int j = 0; j < 4; ++j)
      {
        int   row = std::min(sh - 1, std::max(0, yi - 1 + j));
        float r = 0.f;
        for (int i = 0; i < 4; ++i)
        {
          int col = std::min(sw - 1, std::max(0, xi - 1 + i));
          r += wx[i] * src(col, row);
        }
        val += wy[j] * r;
      }
      dst(x, y) = val;
    }
  }

  return dst;
}

} // namespace hmap
