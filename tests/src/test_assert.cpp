#include <cmath>
#include <limits>

#include "highmap/internal/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(AssertAlmostEqual, RejectsNonfiniteOperandsInEitherOrder)
{
  const Array reference(glm::ivec2(2, 2), 1.f);
  for (float value : {std::numeric_limits<float>::quiet_NaN(),
                      std::numeric_limits<float>::infinity(),
                      -std::numeric_limits<float>::infinity()})
  {
    Array invalid = reference;
    invalid(1, 1) = value;
    AssertResults result;
    EXPECT_FALSE(assert_almost_equal(invalid, reference, 1e-6f, "", &result));
    EXPECT_FALSE(result.ret);
    EXPECT_EQ(result.msg, "non-finite input value");
    EXPECT_FALSE(assert_almost_equal(reference, invalid));
    EXPECT_FALSE(assert_almost_equal(invalid, invalid));
  }
}

TEST(AssertAlmostEqual, RejectsNonfiniteTolerance)
{
  const Array input(glm::ivec2(2, 2), 1.f);
  EXPECT_FALSE(assert_almost_equal(input,
                                   input,
                                   std::numeric_limits<float>::quiet_NaN()));
  EXPECT_FALSE(assert_almost_equal(input,
                                   input,
                                   std::numeric_limits<float>::infinity()));
}

TEST(AssertAlmostEqual, ZeroFieldComparisonIsSymmetric)
{
  const Array zero(glm::ivec2(2, 2));
  const Array one(glm::ivec2(2, 2), 1.f);
  EXPECT_FALSE(assert_almost_equal(zero, one));
  EXPECT_FALSE(assert_almost_equal(one, zero));
  EXPECT_TRUE(assert_almost_equal(zero, zero));
  EXPECT_TRUE(assert_almost_equal(zero, one, 4.f));
  EXPECT_TRUE(assert_almost_equal(one, zero, 4.f));
}

TEST(AssertAlmostEqual, PreservesRelativeL1ErrorAndTolerance)
{
  const Array   a = Array({{1.f, -3.f}, {0.f, 0.f}});
  const Array   b = Array({{1.f, -3.004f}, {0.f, 0.f}});
  AssertResults result;
  EXPECT_FALSE(assert_almost_equal(a, b, 0.0005f, "", &result));
  EXPECT_NEAR(result.diff, 0.001f, 1e-7f);
  EXPECT_FLOAT_EQ(result.count, 0.25f);
  EXPECT_TRUE(assert_almost_equal(a, b, 0.002f));
  EXPECT_TRUE(assert_almost_equal(b, a, 0.002f));
}

TEST(AssertAlmostEqual, FiniteExtremesDoNotOverflowNormalization)
{
  const float max = std::numeric_limits<float>::max();
  const Array large(glm::ivec2(4, 4), max);
  const Array half(glm::ivec2(4, 4), max / 2.f);
  const Array negative(glm::ivec2(4, 4), -max);
  EXPECT_TRUE(assert_almost_equal(large, large));
  EXPECT_FALSE(assert_almost_equal(large, half));
  EXPECT_FALSE(assert_almost_equal(half, large));
  EXPECT_FALSE(assert_almost_equal(large, negative));
  AssertResults result;
  EXPECT_FALSE(assert_almost_equal(large, half, 1e-6f, "", &result));
  EXPECT_FLOAT_EQ(result.diff, 1.f);
}

TEST(AssertAlmostEqual, SuccessfulComparisonResetsReusedResults)
{
  const Array   zero(glm::ivec2(2, 2));
  const Array   one(glm::ivec2(2, 2), 1.f);
  AssertResults result;
  EXPECT_FALSE(assert_almost_equal(zero, one, 1e-6f, "", &result));
  EXPECT_FALSE(result.ret);
  EXPECT_TRUE(assert_almost_equal(one, one, 1e-6f, "", &result));
  EXPECT_TRUE(result.ret);
  EXPECT_TRUE(result.msg.empty());
  EXPECT_FLOAT_EQ(result.diff, 0.f);
  EXPECT_FLOAT_EQ(result.count, 0.f);
}
