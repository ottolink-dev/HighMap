/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU
 * General Public License. The full license is in the LICENSE file. */

#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/range.hpp"

using namespace hmap;

// With FMA contraction (Apple Clang on arm64) 'x * scale + offset' left the
// minimum slightly off vmin, and the maximum can round past vmax on any
// platform. pow-based filters then returned NaN, which spread into black
// tiles in Hesiod on macOS.
TEST(Remap, StaysWithinTargetRange)
{
  std::mt19937                          gen(3);
  std::uniform_real_distribution<float> dis(-1.f, 1.f);

  for (int rep = 0; rep < 200; ++rep)
  {
    Array       a(glm::ivec2(64, 48));
    const float shift = 1000.f * dis(gen);
    const float scale = std::pow(10.f, 3.f * dis(gen));
    for (float &v : a.vector)
      v = shift + scale * dis(gen);

    remap(a, 0.f, 1.f);
    EXPECT_EQ(a.min(), 0.f);
    EXPECT_LE(a.max(), 1.f);
    for (float v : a.vector)
      ASSERT_TRUE(std::isfinite(std::pow(v, 0.8f)));
  }
}

TEST(Remap, ExplicitRangeMapsMinimumExactly)
{
  Array a(glm::ivec2(4, 4));
  for (size_t k = 0; k < a.vector.size(); ++k)
    a.vector[k] = 0.1f + 0.37f * static_cast<float>(k);

  const float from_min = a.vector.front();
  const float from_max = a.vector.back();

  remap(a, 0.f, 1.f, from_min, from_max);
  EXPECT_EQ(a.vector.front(), 0.f);
}

TEST(Gain, FiniteJustOutsideUnitRange)
{
  Array a(glm::ivec2(2, 2));
  a.vector = {-1e-10f, 0.f, 1.f, 1.00000012f};

  gain(a, 1.5f);
  for (float v : a.vector)
  {
    EXPECT_TRUE(std::isfinite(v));
    EXPECT_GE(v, 0.f);
    EXPECT_LE(v, 1.f);
  }
}
