/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/export/image.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"

namespace hmap
{

bool assert_almost_equal(const Array       &a,
                         const Array       &b,
                         float              tolerance,
                         const std::string &fname,
                         AssertResults     *p_results)
{
  // --- check array shapes

  if (a.shape != b.shape)
  {
    std::string msg = "input arrays do not have "
                      "the same shape: a{" +
                      std::to_string(a.shape.x) + ", " +
                      std::to_string(a.shape.y) + "} vs. b {" +
                      std::to_string(b.shape.x) + ", " +
                      std::to_string(b.shape.y) + "}";

    if (p_results)
    {
      p_results->msg = msg;
      p_results->ret = false;
    }

    return false;
  }

  // --- compare

  float diff = 0.f;
  float sum_a = 0.f;
  float sum_b = 0.f;
  int   diff_count = 0;

  {
    for (int j = 0; j < a.shape.y; ++j)
      for (int i = 0; i < a.shape.x; ++i)
      {
        float delta = std::abs(a(i, j) - b(i, j));

        diff += delta;
        sum_a += std::abs(a(i, j));
        sum_b += std::abs(b(i, j));

        if (delta > tolerance) diff_count++;
      }

    float sum = std::min(sum_a, sum_b);

    // if we're comparing zeroes
    if (sum == 0.f)
      diff = sum_a - sum_b;
    else
      diff /= sum;
  }

  if (p_results)
  {
    p_results->diff = diff;
    p_results->tolerance = tolerance;
    p_results->count = (float)diff_count / (float)(a.shape.x * a.shape.y);
  }

  if (fname != "")
  {
    Array ar = a;
    Array br = b;
    Array d = abs(a - b);

    remap(ar);
    remap(br, 0.f, 1.f, a.min(), a.max());
    remap(d);

    export_banner_png(fname, {ar, br, d}, hmap::Cmap::JET);
  }

  if (diff > tolerance)
  {
    std::string msg = "difference above tolerance";

    if (p_results)
    {
      p_results->msg = msg;
      p_results->ret = false;
    }

    return false;
  }

  return true;
}

bool is_subset(const Array &a, const Array &b)
{
  // true if all non-zero pixels of a are also non-zero in b
  for (int j = 0; j < a.shape.y; ++j)
    for (int i = 0; i < a.shape.x; ++i)
      if (a(i, j) != 0.f && b(i, j) == 0.f) return false;
  return true;
}

} // namespace hmap
