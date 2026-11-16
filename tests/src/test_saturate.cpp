#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Saturate, ClampOnly_KZero)
{
  Array input = Array({{0.f, 0.5f, 1.f}});

  Array out = input;
  saturate(out, 0.2f, 0.8f, 0.f, 1.f, 0.f);

  EXPECT_NEAR(out(0, 0), 0.f, 1e-5);
  EXPECT_NEAR(out(0, 1), 0.5f, 1e-5);
  EXPECT_NEAR(out(0, 2), 1.f, 1e-5);
}

TEST(Saturate, RemapRange)
{
  Array input = Array({{0.f, 0.5f, 1.f}});

  Array out = input;
  saturate(out, 10.f, 20.f, 0.f, 1.f, 0.f);

  EXPECT_NEAR(out(0, 0), 0.f, 1e-5);
  EXPECT_NEAR(out(0, 1), 0.f, 1e-5);
  EXPECT_NEAR(out(0, 2), 0.f, 1e-5);
}

TEST(Saturate, MonotonicityPreserved)
{
  Array input = Array({{0.f, 0.2f, 0.4f, 0.6f, 0.8f, 1.f}});

  Array out = input;
  saturate(out, 0.f, 1.f, 0.f, 1.f, 0.f);

  for (int i = 1; i < input.shape.x; ++i)
    EXPECT_LE(out(i - 1, 0), out(i, 0));
}

TEST(Saturate, AutoRangeBasic)
{
  Array input = Array({{0.5f, 2.f, 3.f}});

  Array out = input;
  saturate(out, 0.f, 1.f, 0.f);

  EXPECT_NEAR(out(0, 0), 0.5f, 1e-5);
  EXPECT_NEAR(out(0, 1), 3.f, 1e-5);
  EXPECT_NEAR(out(0, 2), 3.f, 1e-5);
}

TEST(Saturate, AutoRangeClamp)
{
  Array input = Array({{0.f, 5.f, 10.f}});

  Array out = input;
  saturate(out, 2.f, 8.f, 0.f);

  EXPECT_FLOAT_EQ(out(0, 0), 0.f);
  EXPECT_FLOAT_EQ(out(0, 1), 5.f);
  EXPECT_FLOAT_EQ(out(0, 2), 10.f);
}

TEST(Saturate, AutoRangeConstantField)
{
  Array input = Array({{5.f, 5.f, 5.f}});

  Array out = input;
  saturate(out, 0.f, 1.f, 0.f);

  // constant input → constant output
  EXPECT_TRUE(assert_almost_equal(out, Array({{5.f, 5.f, 5.f}})));
}

TEST(SaturatePercentile, BasicPercentile)
{
  Array input = Array({{0.f, 1.f, 2.f, 3.f, 4.f}});

  Array out = input;
  saturate_percentile(out, 0.2f, 0.8f, 0.f);

  // result should be clamped inside percentile range
  float min_val = out.min();
  float max_val = out.max();

  EXPECT_LE(min_val, out(0, 1)); // ~20th percentile
  EXPECT_GE(max_val, out(0, 3)); // ~80th percentile
}

TEST(SaturatePercentile, Monotonicity)
{
  Array input = Array({{0.f, 1.f, 2.f, 3.f, 4.f}});

  Array out = input;
  saturate_percentile(out, 0.1f, 0.9f, 0.f);

  for (int i = 1; i < input.shape.x; ++i)
    EXPECT_LE(out(i - 1, 0), out(i, 0));
}

TEST(Saturate, SmoothClampRangePreserved)
{
  Array input = Array({{-1.f, 0.f, 1.f, 2.f}});

  Array out = input;
  saturate(out, 0.f, 1.f, 0.f, 1.f, 0.5f);

  for (int i = 0; i < out.shape.x; ++i)
  {
    EXPECT_GE(out(i, 0), 0.f);
    EXPECT_LE(out(i, 0), 1.f);
  }
}

TEST(Saturate, SmoothClampMonotonic)
{
  Array input = Array({{0.f, 0.2f, 0.4f, 0.6f, 0.8f, 1.f}});

  Array out = input;
  saturate(out, 0.f, 1.f, 0.f, 1.f, 0.5f);

  for (int i = 1; i < out.shape.x; ++i)
    EXPECT_LE(out(i - 1, 0), out(i, 0));
}
