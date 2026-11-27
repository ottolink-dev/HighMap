/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <cmath>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/export/image.hpp"
#include "highmap/internal/assert.hpp"
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
  // HighMap release builds use fast-math, which can fold std::isfinite to
  // true. Inspect the IEEE-754 exponent bits instead of FP comparisons.
  static_assert(sizeof(float) == sizeof(std::uint32_t) &&
                std::numeric_limits<float>::is_iec559);
  const auto finite = [](float value)
  {
    return (std::bit_cast<std::uint32_t>(value) & 0x7f800000u) != 0x7f800000u;
  };
  const auto invalid = [&](const char *message)
  {
    if (p_results)
    {
      p_results->msg = message;
      p_results->ret = false;
      p_results->diff = std::numeric_limits<float>::infinity();
      p_results->tolerance = tolerance;
    }
    return false;
  };
  if (!finite(tolerance)) return invalid("non-finite tolerance");

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

  // Accumulate in double: even finite float operands can overflow a float
  // difference or sum and turn the final normalized error into NaN.
  double diff = 0.;
  double sum_a = 0.;
  double sum_b = 0.;
  int   diff_count = 0;

  {
    for (int j = 0; j < a.shape.y; ++j)
      for (int i = 0; i < a.shape.x; ++i)
      {
        if (!finite(a(i, j)) || !finite(b(i, j)))
          return invalid("non-finite input value");
        double delta = std::abs((double)a(i, j) - (double)b(i, j));

        diff += delta;
        sum_a += std::abs((double)a(i, j));
        sum_b += std::abs((double)b(i, j));

        if (delta > tolerance) diff_count++;
      }

    double sum = std::min(sum_a, sum_b);

    // if we're comparing zeroes
    if (sum == 0.f)
      diff = std::abs(sum_a - sum_b);
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

  if (p_results)
  {
    p_results->ret = true;
    p_results->msg.clear();
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
