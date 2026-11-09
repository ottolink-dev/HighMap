#include "highmap/array.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(RangeChop, BasicThreshold)
{
  Array a = Array({{0.2f, 0.6f, 1.f}});
  chop(a, 0.5f);

  EXPECT_FLOAT_EQ(a(0, 0), 0.f);
  EXPECT_FLOAT_EQ(a(0, 1), 0.6f);
  EXPECT_FLOAT_EQ(a(0, 2), 1.f);
}

TEST(RangeChop, AllBelowThreshold)
{
  Array a = Array({{0.1f, 0.2f}});
  chop(a, 0.5f);

  EXPECT_EQ(a.min(), 0.f);
}

TEST(RangeChopMaxSmooth, BehaviorRegions)
{
  Array a = Array({{0.f, 0.6f, 1.f}});

  chop_max_smooth(a, 0.5f);

  EXPECT_EQ(a(0, 0), 0.f); // below threshold unchanged-ish
  EXPECT_GE(a(0, 1), 0.f); // mid region modified
}

TEST(RangeClamp, Basic)
{
  Array a = Array({{-1.f, 0.5f, 2.f}});
  clamp(a, 0.f, 1.f);

  EXPECT_FLOAT_EQ(a(0, 0), 0.f);
  EXPECT_FLOAT_EQ(a(0, 1), 0.5f);
  EXPECT_FLOAT_EQ(a(0, 2), 1.f);
}

TEST(RangeClamp, SymmetricMode)
{
  Array a = Array({{-2.f, 0.f, 2.f}});
  clamp(a, 1.f, ClampMode::BOTH);

  EXPECT_LE(a(0, 0), 1.f);
  EXPECT_GE(a(0, 2), -1.f);
}

TEST(RangeClamp, PositiveOnly)
{
  Array a = Array({{-1.f, 2.f}});
  clamp(a, 1.f, ClampMode::POSITIVE_ONLY);

  EXPECT_EQ(a(0, 0), 0.f);
}

TEST(RangeClampMax, Scalar)
{
  Array a = Array({{1.f, 5.f, 10.f}});
  clamp_max(a, 5.f);

  EXPECT_EQ(a(0, 2), 5.f);
}

TEST(RangeClampMin, Scalar)
{
  Array a = Array({{1.f, 5.f, 10.f}});
  clamp_min(a, 5.f);

  EXPECT_EQ(a(0, 0), 5.f);
}

TEST(RangeClampSmooth, NoCrashes)
{
  Array a = Array({{0.f, 0.5f, 1.f}});
  clamp_smooth(a, 0.2f, 0.8f, 0.1f);

  EXPECT_TRUE(std::isfinite(a(0, 1)));
}

TEST(RangeMaximum, ArrayArray)
{
  Array a = Array({{1.f, 5.f}});
  Array b = Array({{2.f, 3.f}});

  Array c = maximum(a, b);

  EXPECT_EQ(c(0, 0), 2.f);
  EXPECT_EQ(c(0, 1), 5.f);
}

TEST(RangeMinimum, ArrayArray)
{
  Array a = Array({{1.f, 5.f}});
  Array b = Array({{2.f, 3.f}});

  Array c = minimum(a, b);

  EXPECT_EQ(c(0, 0), 1.f);
  EXPECT_EQ(c(0, 1), 3.f);
}

TEST(RangeMaxMinSmooth, Consistency)
{
  float a = 1.f, b = 2.f;

  float mx = maximum_smooth(a, b, 0.1f);
  float mn = minimum_smooth(a, b, 0.1f);

  EXPECT_GE(mx, std::max(a, b));
  EXPECT_LE(mn, std::min(a, b));
}

TEST(RangeRemap, NormalizeRange)
{
  Array a = Array({{0.f, 1.f, 2.f}});

  remap(a, 0.f, 1.f);

  EXPECT_NEAR(a.min(), 0.f, 1e-5);
  EXPECT_NEAR(a.max(), 1.f, 1e-5);
}

TEST(RangeRemap, FlatArray)
{
  Array a = Array({{2.f, 2.f, 2.f}});

  remap(a, 0.f, 1.f);

  EXPECT_EQ(a(0, 0), 0.f);
}

TEST(RangeRemap, ExplicitRange)
{
  Array a = Array({{10.f, 20.f, 30.f}});

  remap(a, 0.f, 1.f, 10.f, 30.f);

  EXPECT_NEAR(a(0, 0), 0.f, 1e-5);
  EXPECT_NEAR(a(0, 2), 1.f, 1e-5);
}

TEST(RangeRescale, ScalingOnly)
{
  Array a = Array({{1.f, 2.f, 3.f}});
  rescale(a, 2.f, 0.f);

  EXPECT_EQ(a(0, 0), 2.f);
}

TEST(RangeRescale, WithReference)
{
  Array a = Array({{1.f, 2.f, 3.f}});
  rescale(a, 2.f, 2.f);

  EXPECT_NEAR(a(0, 1), 2.f, 1e-5);
}

TEST(RangeNormalizeAndClip, WithoutSmoothing)
{
  Array a = Array({{-10.f, 0.f, 10.f}});
  normalize_and_clip(a, 1.0f);

  EXPECT_FLOAT_EQ(a(0, 0), -1.f);
  EXPECT_FLOAT_EQ(a(0, 1), 0.f);
  EXPECT_FLOAT_EQ(a(0, 2), 1.f);
}

TEST(RangeNormalizeAndClip, WithSmoothing)
{
  Array a_smooth = Array({{-100.f, -2.f, 0.f, 2.f, 100.f}});
  Array a_raw = a_smooth;

  normalize_and_clip(a_raw, 0.5f, false);
  normalize_and_clip(a_smooth, 0.5f, true);

  // Mean is 0. Both should keep mean at 0
  EXPECT_FLOAT_EQ(a_smooth(0, 2), 0.f);

  // Extremes were clamped to -1 and 1, so smoothstep maps them to -1 and 1
  EXPECT_FLOAT_EQ(a_smooth(0, 0), -1.f);
  EXPECT_FLOAT_EQ(a_smooth(0, 4), 1.f);

  // Smoothstep S-curve softens transitions near boundaries
  // For intermediate positive value, t in (0.5, 1.0), so smoothstep3(t) > t and
  // thus a_smooth > a_raw
  EXPECT_GT(a_smooth(0, 3), a_raw(0, 3));
  EXPECT_LT(a_smooth(0, 1), a_raw(0, 1));

  for (int j = 0; j < 5; ++j)
  {
    EXPECT_GE(a_smooth(0, j), -1.f);
    EXPECT_LE(a_smooth(0, j), 1.f);
  }
}
