/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/convolve.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/kernels.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/random.hpp"
#include "highmap/range.hpp"
#include "highmap/transform.hpp"
#include "highmap/vectors.hpp"

#define NSIGMA 2

namespace hmap
{

void adaptive_relief(Array &array, float strength, float clamp_ratio)
{
  const glm::ivec2 &shape = array.shape;
  Array             result = array;

#pragma omp parallel for collapse(2)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      // 3x3 neighborhood bounds
      int im1 = std::max(0, i - 1);
      int ip1 = std::min(shape.x - 1, i + 1);
      int jm1 = std::max(0, j - 1);
      int jp1 = std::min(shape.y - 1, j + 1);

      float c = array(i, j);
      float n = array(i, jm1);
      float s = array(i, jp1);
      float w = array(im1, j);
      float e = array(ip1, j);
      float nw = array(im1, jm1);
      float ne = array(ip1, jm1);
      float sw = array(im1, jp1);
      float se = array(ip1, jp1);

      // Local min and max
      float min_z = std::min({c, n, s, w, e, nw, ne, sw, se});
      float max_z = std::max({c, n, s, w, e, nw, ne, sw, se});

      // Discrete Laplacian (high-pass filter)
      // Standard 3x3 kernel: [[1, 2, 1], [2, -12, 2], [1, 2, 1]] / 16
      // or [[0.5, 1, 0.5], [1, -6, 1], [0.5, 1, 0.5]]
      float laplacian = c -
                        (0.5f * (n + s + w + e) + 0.25f * (nw + ne + sw + se)) /
                            3.0f;

      // Local contrast / activity across neighborhood
      float contrast = max_z - min_z;

      // Compute edge/slope magnitude to reduce sharpening near extreme cliffs
      float gx = (e - w) * 0.5f;
      float gy = (s - n) * 0.5f;
      float grad_sq = gx * gx + gy * gy;

      // Edge dampening factor (protects steep cliffs/discontinuities)
      float edge_factor = 1.0f / (1.0f + 4.0f * grad_sq /
                                             (contrast * contrast + 1e-7f));

      // Adaptive weight
      float w_adaptive = strength * edge_factor;

      // Sharpened value
      float val = c + w_adaptive * laplacian;

      // Constrain / anti-overshoot against neighborhood extrema
      if (clamp_ratio > 0.f)
      {
        float clamped_val = std::clamp(val, min_z, max_z);
        val = std::lerp(val, clamped_val, clamp_ratio);
      }

      result(i, j) = val;
    }

  array = std::move(result);
}

void adaptive_relief(Array       &array,
                     const Array *p_mask,
                     float        strength,
                     float        clamp_ratio)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { adaptive_relief(a, strength, clamp_ratio); });
}

void equalize(Array &array)
{
  Array flat_ref = hmap::white(array.shape, 0.f, 1.f, 0);
  match_histogram(array, flat_ref);
}

void equalize(Array &array, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [](Array &a) { equalize(a); });
}

void expand(Array &array, int ir, int iterations)
{
  Array array_new = array;
  int   ni = array.shape.x;
  int   nj = array.shape.y;
  Array k = cubic_pulse({2 * ir + 1, 2 * ir + 1});

  for (int it = 0; it < iterations; ++it)
  {
    for (int j = 0; j < nj; j++)
    {
      int q1 = std::max(0, j - ir) - j;
      int q2 = std::min(nj, j + ir + 1) - j;
      for (int i = 0; i < ni; i++)
      {
        int p1 = std::max(0, i - ir) - i;
        int p2 = std::min(ni, i + ir + 1) - i;

        for (int q = q1; q < q2; q++)
          for (int p = p1; p < p2; p++)
          {
            float v = array(i + p, j + q) * k(p + ir, q + ir);
            array_new(i, j) = std::max(array_new(i, j), v);
          }
      }
    }

    array = array_new;
  }
}

void expand(Array &array, int ir, const Array *p_mask, int iterations)
{
  apply_with_mask(array, p_mask, [&](Array &a) { expand(a, ir, iterations); });
}

void expand(Array &array, const Array &kernel, int iterations)
{
  Array array_new = array;
  int   ni = array.shape.x;
  int   nj = array.shape.y;

  int ri1 = (int)(0.5f * kernel.shape.x);
  int ri2 = kernel.shape.x - ri1 - 1;
  int rj1 = (int)(0.5f * kernel.shape.y);
  int rj2 = kernel.shape.y - rj1 - 1;

  for (int it = 0; it < iterations; ++it)
  {
    for (int j = 0; j < nj; j++)
    {
      int q1 = std::max(0, j - rj1) - j;
      int q2 = std::min(nj, j + rj2 + 1) - j;
      for (int i = 0; i < ni; i++)
      {
        int p1 = std::max(0, i - ri1) - i;
        int p2 = std::min(ni, i + ri2 + 1) - i;

        for (int q = q1; q < q2; q++)
          for (int p = p1; p < p2; p++)
          {
            float v = array(i + p, j + q) * kernel(p + ri1, q + rj1);
            array_new(i, j) = std::max(array_new(i, j), v);
          }
      }
    }
    array = array_new;
  }
}

void expand(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { expand(a, kernel, iterations); });
}

void expand_directional(Array       &array,
                        int          ir,
                        float        angle,
                        float        aspect_ratio,
                        float        anisotropy,
                        const Array *p_mask)
{
  Array kernel = cubic_pulse_directional(glm::ivec2(2 * ir + 1, 2 * ir + 1),
                                         angle,
                                         aspect_ratio,
                                         anisotropy);
  expand(array, kernel, p_mask);
}

void fold(Array &array, int iterations, float k)
{
  fold(array, array.min(), array.max(), iterations, k);
}

void fold(Array &array, float vmin, float vmax, int iterations, float k)
{
  array -= vmin;
  float vref = (vmax - vmin) / (iterations + 1.f);

  for (int it = 0; it < iterations; it++)
  {
    array -= vref;

    if (k == 0.f)
      array = abs(array);
    else
      array = abs_smooth(array, k);
  }
}

void gain(Array &array, float factor)
{
  auto lambda = [&factor](float x)
  {
    return x < 0.5 ? 0.5f * std::pow(2.f * x, factor)
                   : 1.f - 0.5f * std::pow(2.f * (1.f - x), factor);
  };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void gain(Array &array, float factor, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { gain(a, factor); });
}

void gamma_correction(Array &array, float gamma)
{
  auto lambda = [&gamma](float x) { return std::pow(x, gamma); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void gamma_correction(Array &array, float gamma, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { gamma_correction(a, gamma); });
}

void gamma_correction_local(Array &array, float gamma, int ir, float k)
{
  Array amin = local_min(array, ir);
  Array amax = local_max(array, ir);

  smooth_cpulse(amin, ir);
  smooth_cpulse(amax, ir);

  if (k != 0) // with smoothing
  {
    for (int j = 0; j < array.shape.y; j++)
      for (int i = 0; i < array.shape.x; i++)
      {
        float v = std::abs(array(i, j) - amin(i, j)) /
                  (amax(i, j) - amin(i, j) + 1e-30);
        v = std::sqrt(v * v + k);
        array(i, j) = std::pow(v, gamma) * (amax(i, j) - amin(i, j)) +
                      amin(i, j);
      }
  }
  else // without smoothing
  {
    for (int j = 0; j < array.shape.y; j++)
      for (int i = 0; i < array.shape.x; i++)
      {
        float v = std::abs(array(i, j) - amin(i, j)) /
                  (amax(i, j) - amin(i, j) + 1e-30);
        array(i, j) = std::pow(v, gamma) * (amax(i, j) - amin(i, j)) +
                      amin(i, j);
      }
  }
}

void gamma_correction_local(Array       &array,
                            float        gamma,
                            int          ir,
                            const Array *p_mask,
                            float        k)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { gamma_correction_local(a, gamma, ir, k); });
}

void kuwahara(Array &array, int ir, float mix_ratio)
{

  Array array_buffered = generate_buffered_array(array,
                                                 glm::ivec4(ir, ir, ir, ir));
  Array array_out(array_buffered.shape);

  for (int j = ir; j < array_buffered.shape.y - ir; j++)
    for (int i = ir; i < array_buffered.shape.x - ir; i++)
    {
      // build quadrants
      Array q1 = array_buffered.extract_slice(
          glm::ivec4(i - ir, i + 1, j - ir, j + 1));
      Array q2 = array_buffered.extract_slice(
          glm::ivec4(i - ir, i + 1, j + 1, j + ir));
      Array q3 = array_buffered.extract_slice(
          glm::ivec4(i + 1, i + ir, j - ir, j + 1));
      Array q4 = array_buffered.extract_slice(
          glm::ivec4(i + 1, i + ir, j + 1, j + ir));

      std::vector<float> means = {q1.mean(), q2.mean(), q3.mean(), q4.mean()};
      std::vector<float> stds = {q1.std(), q2.std(), q3.std(), q4.std()};

      size_t imin = std::min_element(stds.begin(), stds.end()) - stds.begin();
      array_out(i, j) = means[imin];
    }

  if (mix_ratio == 1.f)
  {
    array = array_out.extract_slice(glm::ivec4(ir,
                                               array_buffered.shape.x - ir,
                                               ir,
                                               array_buffered.shape.y - ir));
  }
  else
  {
    array = lerp(
        array,
        array_out.extract_slice(glm::ivec4(ir,
                                           array_buffered.shape.x - ir,
                                           ir,
                                           array_buffered.shape.y - ir)),
        mix_ratio);
  }
}

void kuwahara(Array &array, int ir, const Array *p_mask, float mix_ratio)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { kuwahara(a, ir, p_mask ? 1.f : mix_ratio); });
}

void laplace(Array &array, float sigma, int iterations)
{
  glm::ivec2 shape = array.shape;

  for (int it = 0; it < iterations; it++)
  {
    Array delta(shape);

    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        float center = array(i, j);
        float sum_neighbors = 0.0f;

        // Accumulate 8 neighbors with boundary clamping
        for (int dj = -1; dj <= 1; dj++)
          for (int di = -1; di <= 1; di++)
          {
            if (di == 0 && dj == 0) continue;

            int ni = std::clamp(i + di, 0, shape.x - 1);
            int nj = std::clamp(j + dj, 0, shape.y - 1);
            sum_neighbors += array(ni, nj);
          }

        // 8-neighbor Laplacian
        delta(i, j) = sum_neighbors - 8.0f * center;
      }

    // Apply diffusion update
    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
        array(i, j) += sigma * delta(i, j);
  }
}

void laplace(Array &array, const Array *p_mask, float sigma, int iterations)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { laplace(a, sigma, iterations); });
}

void laplace_edge_preserving(Array &array,
                             float  talus,
                             float  sigma,
                             int    iterations)
{
  for (int it = 0; it < iterations; it++)
  {
    Array c = gradient_norm(array);
    c = 1.f / (1.f + c * c / (talus * talus));

    Array dcx = gradient_x(c);
    Array dcy = gradient_y(c);
    Array dzx = gradient_x(array);
    Array dzy = gradient_y(array);
    Array delta = laplacian(array);

    array += sigma * (dcx * dzx + dcy * dzy + c * delta);
  }
}

void laplace_edge_preserving(Array       &array,
                             float        talus,
                             const Array *p_mask,
                             float        sigma,
                             int          iterations)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { laplace_edge_preserving(a, talus, sigma, iterations); });
}

void low_pass_high_order(Array &array, int order, float sigma)
{
  Array df = array;

  // filtering coefficients
  std::vector<float> kernel;

  switch (order)
  {
  case (5): kernel = {0.0625f, -0.25f, 0.375f, -0.25f, 0.0625f}; break;

  case (7):
    kernel = {-0.015625f,
              0.09375f,
              -0.234375f,
              0.3125f,
              -0.234375f,
              0.09375f,
              -0.015625f};
    break;

  case (9):
    kernel = {0.00390625f,
              -0.03125f,
              0.109375f,
              -0.21875f,
              0.2734375f,
              -0.21875f,
              0.109375f,
              -0.03125f,
              0.00390625f};
    break;
  }

  df = convolve1d_i(array, kernel);
  array -= sigma * df;

  df = convolve1d_j(array, kernel);
  array -= sigma * df;
}

void make_binary(Array &array, float threshold)
{
  auto lambda = [&threshold](float a) { return a > threshold ? 1.f : 0.f; };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 array.vector.begin(),
                 lambda);
}

void match_histogram(Array &array, const Array &array_reference)
{
  std::vector<size_t> ki = argsort(array.vector);
  std::vector<size_t> kr = argsort(array_reference.vector);

  size_t n = ki.size();
  size_t nr = kr.size();

  for (size_t i = 0; i < n; i++)
  {
    // map rank i in source → fractional rank in reference
    float  t = static_cast<float>(i) / (n - 1) * (nr - 1);
    size_t j0 = static_cast<size_t>(t);
    size_t j1 = std::min(j0 + 1, nr - 1);
    float  f = t - static_cast<float>(j0);

    array.vector[ki[i]] = (1.f - f) * array_reference.vector[kr[j0]] +
                          f * array_reference.vector[kr[j1]];
  }
}

Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 int          iterations,
                 bool         talus_weighted)
{
  const glm::ivec2 shape = array.shape;
  Array            array_next = Array(shape);
  Array            array_prev = array;

  for (int it = 0; it < iterations; it++)
  {

    if (talus_weighted)
    {
      for (int j = 0; j < shape.y; j++)
        for (int i = 0; i < shape.x; i++)
        {
          float sum = 0.f;
          float norm = 0.f;

          for (int q = j - ir; q < j + ir + 1; ++q)
            for (int p = i - ir; p < i + ir + 1; ++p)
              if (p > 0 && p < shape.x && q > 0 && q < shape.y)
              {
                float dv = std::abs(array_prev(i, j) - array_prev(p, q));
                if (dv < talus)
                {
                  float weight = 1.f - dv / talus;
                  sum += array(p, q) * weight;
                  norm += weight;
                }
              }

          array_next(i, j) = sum / norm;
        }
    }
    else
    {
      for (int j = 0; j < shape.y; j++)
        for (int i = 0; i < shape.x; i++)
        {
          int   count = 0.f;
          float sum = 0.f;

          for (int q = j - ir; q < j + ir + 1; ++q)
            for (int p = i - ir; p < i + ir + 1; ++p)
              if (p > 0 && p < shape.x && q > 0 && q < shape.y)
              {
                float dv = std::abs(array_prev(i, j) - array_prev(p, q));
                if (dv < talus)
                {
                  sum += array(p, q);
                  count++;
                }
              }

          array_next(i, j) = sum / static_cast<float>(count);
        }
    }

    if (iterations > 1) array_prev = array_next;
  }

  return array_next;
}

Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 const Array *p_mask,
                 int          iterations,
                 bool         talus_weighted)
{
  return transform_with_mask(
      array,
      p_mask,
      [&](const Array &a)
      { return mean_shift(a, ir, talus, iterations, talus_weighted); });
}

void median_3x3(Array &array)
{
  Array array_out = Array(array.shape);

  std::vector<float> v(9);

  for (int j = 1; j < array.shape.y - 1; j++)
    for (int i = 1; i < array.shape.x - 1; i++)
    {
      v[0] = array(i - 1, j - 1);
      v[1] = array(i - 1, j);
      v[2] = array(i - 1, j + 1);
      v[3] = array(i, j - 1);
      v[4] = array(i, j);
      v[5] = array(i, j + 1);
      v[6] = array(i + 1, j - 1);
      v[7] = array(i + 1, j);
      v[8] = array(i + 1, j + 1);

      std::sort(v.begin(), v.end());
      array_out(i, j) = v[4];
    }

  for (int j = 1; j < array.shape.y - 1; j++)
    for (int i = 1; i < array.shape.x - 1; i++)
      array(i, j) = array_out(i, j);
}

void median_3x3(Array &array, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [](Array &a) { median_3x3(a); });
}

Array median_pseudo(const Array &array, int ir)
{
  return (local_min(array, ir) + local_max(array, ir) + local_mean(array, ir)) /
         3.f;
}

void normal_displacement(Array &array, float amount, int ir, bool reverse)
{
  Array array_f = array;
  Array array_new = Array(array.shape);

  if (ir > 0) smooth_cpulse(array_f, ir);

  // add a shape factor to avoid artifacts close to the boundaries
  Array factor = biweight(array.shape); // smooth_cosine(array.shape);

  if (reverse) amount = -amount;

  for (int j = 1; j < array.shape.y - 1; j++)
    for (int i = 1; i < array.shape.x - 1; i++)
    {
      glm::vec3 n = array_f.get_normal_at(i, j);

      float x = (float)i - amount * array.shape.x * n.x * factor(i, j);
      float y = (float)j - amount * array.shape.y * n.y * factor(i, j);

      // bilinear interpolation parameters
      int ip = std::clamp((int)x, 0, array.shape.x - 1);
      int jp = std::clamp((int)y, 0, array.shape.y - 1);

      float u = std::clamp(x - (float)ip, 0.f, 1.f);
      float v = std::clamp(y - (float)jp, 0.f, 1.f);

      array_new(i, j) = array.get_value_bilinear_at(ip, jp, u, v);
    }
  fill_borders(array_new);

  array = array_new;
}

void normal_displacement(Array       &array,
                         const Array *p_mask,
                         float        amount,
                         int          ir,
                         bool         reverse)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { normal_displacement(a, amount, ir, reverse); });
}

void plateau(Array &array, int ir, float factor)
{
  Array amin = local_min(array, ir);
  Array amax = local_max(array, ir);

  smooth_cpulse(amin, ir);
  smooth_cpulse(amax, ir);

  array = (array - amin) / (amax - amin + std::numeric_limits<float>::min());
  clamp(array); // keep things under control...
  gain(array, factor);
  array = amin + (amax - amin) * array;
}

void plateau(Array &array, const Array *p_mask, int ir, float factor)
{
  apply_with_mask(array, p_mask, [&](Array &a) { plateau(a, ir, factor); });
}

void reverse_above_theshold(Array       &array,
                            const Array &threshold,
                            float        scaling,
                            float        transition_extent)
{
  for (int j = 0; j < array.shape.y; ++j)
    for (int i = 0; i < array.shape.x; ++i)
    {
      float base_value = array(i, j);
      float limit_value = threshold(i, j);
      float gap = base_value - limit_value;

      if (gap > 0.f)
      {
        float t = std::clamp(gap / transition_extent, 0.f, 1.f);
        t = smoothstep3(t);

        array(i, j) = limit_value - scaling * t * gap;
      }
      else
      {
        float t = std::clamp(-gap / transition_extent, 0.f, 1.f);
        t = smoothstep3(t);

        array(i, j) = lerp(limit_value, base_value, t);
      }
    }
}

void reverse_above_theshold(Array &array,
                            float  threshold,
                            float  scaling,
                            float  transition_extent)
{
  Array threshold_array = Array(array.shape, threshold);

  reverse_above_theshold(array, threshold_array, scaling, transition_extent);
}

void reverse_above_theshold(Array       &array,
                            const Array &threshold,
                            const Array *p_mask,
                            float        scaling,
                            float        transition_extent)
{
  apply_with_mask(
      array,
      p_mask,
      [&](Array &a)
      { reverse_above_theshold(a, threshold, scaling, transition_extent); });
}

void reverse_above_theshold(Array       &array,
                            float        threshold,
                            const Array *p_mask,
                            float        scaling,
                            float        transition_extent)
{
  apply_with_mask(
      array,
      p_mask,
      [&](Array &a)
      { reverse_above_theshold(a, threshold, scaling, transition_extent); });
}

void sharpen(Array &array, float ratio)
{
  Array lp = Array(array.shape);

  for (int j = 1; j < array.shape.y - 1; j++)
    for (int i = 1; i < array.shape.x - 1; i++)
    {
      lp(i, j) = 5.f * array(i, j) - array(i + 1, j) - array(i - 1, j) -
                 array(i, j - 1) - array(i, j + 1);
    }
  extrapolate_borders(lp);
  array = (1.f - ratio) * array + ratio * lp;
}

void sharpen(Array &array, const Array *p_mask, float ratio)
{
  apply_with_mask(array, p_mask, [&](Array &a) { sharpen(a, ratio); });
}

void sharpen_cone(Array &array, int ir, float scale)
{
  Array array_low_pass = array;
  smooth_cone(array_low_pass, ir);
  array += scale * (array - array_low_pass);
}

void sharpen_cone(Array &array, const Array *p_mask, int ir, float scale)
{
  apply_with_mask(array, p_mask, [&](Array &a) { sharpen_cone(a, ir, scale); });
}

void shrink(Array &array, int ir, int iterations)
{
  float amax = array.max();
  array = amax - array;
  expand(array, ir, iterations);
  array = amax - array;
}

void shrink(Array &array, int ir, const Array *p_mask, int iterations)
{
  apply_with_mask(array, p_mask, [&](Array &a) { shrink(a, ir, iterations); });
}

void shrink(Array &array, const Array &kernel, int iterations)
{
  float amax = array.max();
  array = amax - array;
  expand(array, kernel, iterations);
  array = amax - array;
}

void shrink(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { shrink(a, kernel, iterations); });
}

void shrink_directional(Array       &array,
                        int          ir,
                        float        angle,
                        float        aspect_ratio,
                        float        anisotropy,
                        const Array *p_mask)
{
  Array kernel = cubic_pulse_directional(glm::ivec2(2 * ir + 1, 2 * ir + 1),
                                         angle,
                                         aspect_ratio,
                                         anisotropy);
  shrink(array, kernel, p_mask);
}

void smooth_cone(Array &array, int ir)
{
  // define kernel
  const int          nk = 2 * ir + 1;
  std::vector<float> k(nk);

  float sum = 0.f;
  float x0 = (float)nk / 2.f;
  for (int i = 0; i < nk; i++)
  {
    float x = std::abs((float)i - x0) / (float)ir;
    k[i] = 1.f - x;
    sum += k[i];
  }

  // normalize
  for (int i = 0; i < nk; i++)
  {
    k[i] /= sum;
  }

  // eventually convolve
  array = convolve1d_i(array, k);
  array = convolve1d_j(array, k);
}

void smooth_cone(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { smooth_cone(a, ir); });
}

void smooth_cpulse(Array &array, int ir)
{
  // define kernel
  const int          nk = 2 * ir + 1;
  std::vector<float> k1d(nk);

  float sum = 0.f;
  float x0 = (float)ir;
  for (int i = 0; i < nk; i++)
  {
    float x = std::abs((float)i - x0) / (float)ir;
    k1d[i] = std::exp(-0.5f * x * x * 9.f); // σ ≈ ir/3
    sum += k1d[i];
  }

  // normalize
  for (int i = 0; i < nk; i++)
  {
    k1d[i] /= sum;
  }

  // eventually convolve
  array = convolve1d_i(array, k1d);
  array = convolve1d_j(array, k1d);
}

void smooth_cpulse(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { smooth_cpulse(a, ir); });
}

void smooth_cpulse_edge_removing(Array &array,
                                 float  talus,
                                 float  talus_width,
                                 int    ir)
{
  Array c = gradient_norm(array);
  c = sigmoid(c, talus_width, 0.f /* vmin */, 1.f /* vmax */, talus);
  expand(c, ir);
  smooth_cpulse(array, ir, &c);
}

void smooth_flat(Array &array, int ir)
{
  // define kernel
  const int          nk = 2 * ir + 1;
  std::vector<float> k(nk);

  std::fill(k.begin(), k.end(), 1.f / (2.f * nk + 1.f));

  // eventually convolve
  array = convolve1d_i(array, k);
  array = convolve1d_j(array, k);
}

void smooth_gaussian(Array &array, int ir)
{
  // define Gaussian kernel (we keep NSIGMA standard deviations of the
  // kernel support)
  const int          nk = NSIGMA * (2 * ir + 1);
  std::vector<float> k(nk);

  float sum = 0.f;
  float sig2 = (float)(ir * ir);
  float x0 = (float)nk / 2.f;

  for (int i = 0; i < nk; i++)
  {
    float x = (float)i - x0;
    k[i] = std::exp(-0.5f * std::pow(x, 2.f) / sig2);
    sum += k[i];
  }

  // normalize
  for (int i = 0; i < nk; i++)
  {
    k[i] /= sum;
  }

  // eventually convolve
  array = convolve1d_i(array, k);
  array = convolve1d_j(array, k);
}

void smooth_gaussian(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { smooth_gaussian(a, ir); });
}

void smooth_fill(Array &array, int ir, float k, Array *p_deposition_map)
{
  // keep a backup of the input for the deposition map
  Array array_bckp = Array();
  if (p_deposition_map != nullptr) array_bckp = array;

  // smooth filling
  Array array_smooth = array;
  smooth_cpulse(array_smooth, ir);
  array = maximum_smooth(array, array_smooth, k);

  // update map
  if (p_deposition_map)
  {
    *p_deposition_map = array - array_bckp;
    clamp_min(*p_deposition_map, 0.f);
  }
}

void smooth_fill(Array       &array,
                 int          ir,
                 const Array *p_mask,
                 float        k,
                 Array       *p_deposition_map)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { smooth_fill(a, ir, k, p_deposition_map); });
}

void smooth_fill_holes(Array &array, int ir)
{
  Array array_smooth = array;
  smooth_cpulse(array_smooth, ir);

  // mask based on concave regions
  Array mask = curvature_quadric(array_smooth, 0, CurvatureType::CT_MEAN);
  clamp_min(mask, 0.f);
  make_binary(mask);

  int ic = (int)((float)ir / 2.f);
  if (ic > 1) smooth_cpulse(mask, ic);

  array = lerp(array, array_smooth, mask);
}

void smooth_fill_holes(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { smooth_fill_holes(a, ir); });
}

void smooth_fill_smear_peaks(Array &array, int ir)
{
  Array array_smooth = array;
  smooth_cpulse(array_smooth, ir);

  // mask based on concave regions
  Array mask = -curvature_quadric(array_smooth, 0, CurvatureType::CT_MEAN);
  clamp_min(mask, 0.f);
  make_binary(mask);

  int ic = (int)((float)ir / 2.f);
  if (ic > 0) smooth_cpulse(mask, ic);

  array = lerp(array, array_smooth, mask);
}

void smooth_fill_smear_peaks(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { smooth_fill_smear_peaks(a, ir); });
}

void smoothstep_local(Array &array, int ir)
{
  Array amin = local_min(array, ir);
  Array amax = local_max(array, ir);

  smooth_cpulse(amin, ir);
  smooth_cpulse(amax, ir);

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      float v = (array(i, j) - amin(i, j)) / (amax(i, j) - amin(i, j) + 1e-30);
      array(i, j) = smoothstep3(v) * (amax(i, j) - amin(i, j)) + amin(i, j);
    }
}

void smoothstep_local(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { smoothstep_local(a, ir); });
}

void steepen(Array &array, float scale, int ir)
{
  Array dx = gradient_x(array) * ((float)array.shape.x * -scale);
  Array dy = gradient_y(array) * ((float)array.shape.y * -scale);

  smooth_cpulse(dx, ir);
  smooth_cpulse(dy, ir);

  warp(array, &dx, &dy);
}

void steepen(Array &array, float scale, const Array *p_mask, int ir)
{
  apply_with_mask(array, p_mask, [&](Array &a) { steepen(a, scale, ir); });
}

void steepen_convective(Array &array,
                        float  angle,
                        int    iterations,
                        int    ir,
                        float  dt)
{
  Array dx = Array(array.shape);
  Array dy = Array(array.shape);
  float alpha = angle / 180.f * M_PI;
  float ca = std::cos(alpha);
  float sa = std::sin(alpha);

  for (int it = 0; it < iterations; it++)
  {
    if (ir > 0)
    {
      Array array_filtered = array;
      smooth_cpulse(array_filtered, ir);
      gradient_x(array_filtered, dx);
      gradient_y(array_filtered, dy);
    }
    else
    {
      gradient_x(array, dx);
      gradient_y(array, dy);
    }

    dx *= (float)array.shape.x;
    dy *= (float)array.shape.y;
    array *= 1.f - dt * (ca * dx + sa * dy); // == du / dt = - u * du / dx
  }
}

void steepen_convective(Array       &array,
                        float        angle,
                        const Array *p_mask,
                        int          iterations,
                        int          ir,
                        float        dt)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { steepen_convective(a, angle, iterations, ir, dt); });
}

void terrace(Array        &array,
             std::uint32_t seed,
             int           nlevels,
             float         gain,
             float         noise_ratio,
             const Array  *p_noise,
             float         vmin,
             float         vmax)
{
  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(-noise_ratio, noise_ratio);

  // redefine min/max if sentinels values are detected
  if (vmax < vmin)
  {
    vmin = array.min();
    vmax = array.max();
  }

  // defines levels
  std::vector<float> levels = linspace(vmin, vmax, nlevels + 1);
  float              delta = (vmax - vmin) / (float)nlevels;

  // add noise, except for the first and last levels
  for (size_t k = 1; k < levels.size() - 1; k++)
    levels[k] += dis(gen) * delta;

  // apply a gain like filter
  auto lambda = [gain, &levels, vmin, vmax](float x, float noise = 0.f)
  {
    // find level interval
    float y = x + noise;
    y = std::clamp(y, vmin, vmax);

    size_t n = 1;
    while (y > levels[n] && n < levels.size())
      n++;
    n--;

    // rescale value to [0, 1]
    y = (y - levels[n]) / (levels[n + 1] - levels[n]);

    // apply gain
    float gain_clamped = std::clamp(gain, 0.f, 1.f);

    y -= std::atan(-gain_clamped * std::sin(2.f * M_PI * y) /
                   (1.f - gain_clamped * std::cos(2.f * M_PI * y))) /
         M_PI;

    // rescale back to original ammplitude interval
    return y * (levels[n + 1] - levels[n]) + levels[n] - noise;
  };

  if (p_noise)
    std::transform(array.vector.begin(),
                   array.vector.end(),
                   p_noise->vector.begin(),
                   array.vector.begin(),
                   lambda);
  else
    std::transform(array.vector.begin(),
                   array.vector.end(),
                   array.vector.begin(),
                   lambda);
}

void terrace(Array        &array,
             std::uint32_t seed,
             int           nlevels,
             const Array  *p_mask,
             float         gain,
             float         noise_ratio,
             const Array  *p_noise,
             float         vmin,
             float         vmax)
{
  apply_with_mask(
      array,
      p_mask,
      [&](Array &a)
      { terrace(a, seed, nlevels, gain, noise_ratio, p_noise, vmin, vmax); });
}

} // namespace hmap
