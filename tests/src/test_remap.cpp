#include "highmap/array.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/range.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Remap, AutoRangeBasic)
{
  Array input = Array({{0.f, 5.f, 10.f}});
  Array out = input;
  remap(out, 0.f, 1.f);

  EXPECT_NEAR(out(0, 0), 0.f, 1e-5);
  EXPECT_NEAR(out(1, 0), 0.5f, 1e-5);
  EXPECT_NEAR(out(2, 0), 1.f, 1e-5);
}

TEST(Remap, AutoRangeCustomTarget)
{
  Array input = Array({{0.f, 5.f, 10.f}});
  Array out = input;
  remap(out, 10.f, 20.f);

  EXPECT_NEAR(out(0, 0), 10.f, 1e-5);
  EXPECT_NEAR(out(1, 0), 15.f, 1e-5);
  EXPECT_NEAR(out(2, 0), 20.f, 1e-5);
}

TEST(Remap, AutoRangeNegativeValues)
{
  Array input = Array({{-5.f, 0.f, 5.f}});
  Array out = input;
  remap(out, 0.f, 1.f);

  EXPECT_NEAR(out(0, 0), 0.f, 1e-5);
  EXPECT_NEAR(out(1, 0), 0.5f, 1e-5);
  EXPECT_NEAR(out(2, 0), 1.f, 1e-5);
}

TEST(Remap, AutoRangeConstantArray)
{
  Array input = Array({{3.f, 3.f, 3.f}});
  Array out = input;
  remap(out, 0.f, 1.f);

  // degenerate case → filled with vmin
  EXPECT_TRUE(assert_almost_equal(out, Array({{0.f, 0.f, 0.f}})));
}

TEST(Remap, AutoRangeMonotonicity)
{
  Array input = Array({{1.f, 2.f, 3.f, 4.f}});
  Array out = input;
  remap(out, 0.f, 1.f);

  for (int i = 1; i < input.shape.x; ++i)
    EXPECT_LE(out(i - 1, 0), out(i, 0));
}

TEST(Remap, ExplicitRangeBasic)
{
  Array input = Array({{0.f, 5.f, 10.f}});
  Array out = input;
  remap(out, 0.f, 1.f, 0.f, 10.f);

  EXPECT_NEAR(out(0, 0), 0.f, 1e-5);
  EXPECT_NEAR(out(1, 0), 0.5f, 1e-5);
  EXPECT_NEAR(out(2, 0), 1.f, 1e-5);
}

TEST(Remap, ExplicitRangeCustomTarget)
{
  Array input = Array({{0.f, 5.f, 10.f}});
  Array out = input;
  remap(out, 10.f, 20.f, 0.f, 10.f);

  EXPECT_NEAR(out(0, 0), 10.f, 1e-5);
  EXPECT_NEAR(out(1, 0), 15.f, 1e-5);
  EXPECT_NEAR(out(2, 0), 20.f, 1e-5);
}

TEST(Remap, ExplicitRangeOutsideInput)
{
  Array input = Array({{-5.f, 0.f, 5.f}});
  Array out = input;
  remap(out, 0.f, 1.f, 0.f, 10.f);

  // values outside source range are extrapolated
  EXPECT_LT(out(0, 0), 0.f); // -5 maps below 0
  EXPECT_NEAR(out(1, 0), 0.f, 1e-5);
  EXPECT_NEAR(out(2, 0), 0.5f, 1e-5);
}

TEST(Remap, ExplicitRangeConstantSource)
{
  Array input = Array({{5.f, 5.f, 5.f}});
  Array out = input;
  remap(out, 0.f, 1.f, 5.f, 5.f);

  // degenerate case → filled with vmin
  EXPECT_TRUE(assert_almost_equal(out, Array({{0.f, 0.f, 0.f}})));
}

TEST(Remap, ExplicitRangeMonotonicity)
{
  Array input = Array({{1.f, 2.f, 3.f, 4.f}});
  Array out = input;
  remap(out, 0.f, 1.f, 0.f, 5.f);

  for (int i = 1; i < input.shape.x; ++i)
    EXPECT_LE(out(i - 1, 0), out(i, 0));
}

TEST(Remap, RangeMatchesTarget)
{
  Array input = Array({{1.f, 2.f, 3.f, 4.f}});
  Array out = input;
  remap(out, -1.f, 1.f);

  EXPECT_NEAR(out.min(), -1.f, 1e-5);
  EXPECT_NEAR(out.max(), 1.f, 1e-5);
}

TEST(Remap, IdentityMapping)
{
  Array input = Array({{0.f, 0.5f, 1.f}});
  Array out = input;
  remap(out, 0.f, 1.f);

  EXPECT_TRUE(assert_almost_equal(out, input));
}
