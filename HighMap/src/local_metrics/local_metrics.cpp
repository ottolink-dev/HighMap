/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/convolve.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"

namespace hmap
{

Array local_max(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  const int nx = array.shape.x;
  const int ny = array.shape.y;

  Array array_tmp = Array(array.shape);
  Array array_out = Array(array.shape);

  // --- Pass 1: centered sliding max along rows

#pragma omp parallel for schedule(static)
  for (int j = 0; j < ny; ++j)
  {
    std::deque<int> dq;

    // iterate beyond nx-1 to flush the right side of centered windows
    for (int i = 0; i < nx + ir; ++i)
    {
      // add element i if within bounds
      if (i < nx)
      {
        while (!dq.empty() && array(dq.back(), j) <= array(i, j))
          dq.pop_back();
        dq.push_back(i);
      }

      // window for output p = i-ir is [p-ir, p+ir] = [i-2*ir, i]
      while (!dq.empty() && dq.front() < i - 2 * ir)
        dq.pop_front();

      // write output with ir-step delay
      const int p = i - ir;
      if (p >= 0 && p < nx && !dq.empty())
        array_tmp(p, j) = array(dq.front(), j);
    }
  }

  // --- Pass 2: centered sliding max along columns

#pragma omp parallel for schedule(static)
  for (int i = 0; i < nx; ++i)
  {
    std::deque<int> dq;

    for (int j = 0; j < ny + ir; ++j)
    {
      if (j < ny)
      {
        while (!dq.empty() && array_tmp(i, dq.back()) <= array_tmp(i, j))
          dq.pop_back();
        dq.push_back(j);
      }

      while (!dq.empty() && dq.front() < j - 2 * ir)
        dq.pop_front();

      const int p = j - ir;
      if (p >= 0 && p < ny && !dq.empty())
        array_out(i, p) = array_tmp(i, dq.front());
    }
  }

  return array_out;
}

Array local_mean(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = Array(array.shape);

  std::vector<float> k1d(2 * ir + 1);
  for (auto &v : k1d)
    v = 1.f / (float)(2 * ir + 1);

  array_out = convolve1d_i(array, k1d);
  array_out = convolve1d_j(array_out, k1d);

  return array_out;
}

Array local_min(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  return -local_max(-array, ir);
}

Array local_median_deviation(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array mean = local_mean(array, ir);
  Array med = median_pseudo(array, ir); // TODO exact
  return abs(mean - med);
}

Array relative_elevation(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array amin = local_min(array, ir);
  Array amax = local_max(array, ir);

  return (array - amin) / (amax - amin + std::numeric_limits<float>::min());
}

Array roughness(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array smoothed = array;
  smooth_cpulse(smoothed, ir);
  return abs(array - smoothed);
}

Array ruggedness(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array rg(array.shape);

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      int i1 = std::max(i - ir, 0);
      int i2 = std::min(i + ir + 1, array.shape.x);
      int j1 = std::max(j - ir, 0);
      int j2 = std::min(j + ir + 1, array.shape.y);

      for (int p = i1; p < i2; p++)
        for (int q = j1; q < j2; q++)
        {
          float delta = array(i, j) - array(p, q);
          rg(i, j) += delta * delta;
        }
    }

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      rg(i, j) = std::sqrt(rg(i, j));

  return rg;
}

Array rugosity(const Array &z, int ir, bool convex)
{
  if (!validate_non_empty(z)) return Array();

  hmap::Array z_avg = Array(z.shape);
  hmap::Array z_std = Array(z.shape);
  hmap::Array z_skw = Array(z.shape);

  // pre high-pass filter to remove low wavenumbers
  Array zf = z;
  smooth_cpulse(zf, 2 * ir);
  zf = z - zf;

  // use Gaussian windowing instead of a real arithmetic averaging to
  // limit boundary artifacts
  z_avg = zf;
  smooth_cpulse(z_avg, ir);

  z_std = (zf - z_avg) * (zf - z_avg);
  smooth_cpulse(z_std, ir);

  // Fisher-Pearson coefficient of skewness
  z_skw = (zf - z_avg) * (zf - z_avg) * (zf - z_avg);

  // do not filter, surprisingly yields much better results
  // smooth_cpulse(z_skw, ir);

  float tol = 1e-30f * z.ptp();

  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
      if (z_std(i, j) > tol)
        z_skw(i, j) /= std::pow(z_std(i, j), 1.5f);
      else
        z_skw(i, j) = 0.f;

  // keep only "bumpy" rugosities
  if (convex)
    clamp_min(z_skw, 0.f);
  else
    clamp_max(z_skw, 0.f);

  return z_skw;
}

} // namespace hmap
