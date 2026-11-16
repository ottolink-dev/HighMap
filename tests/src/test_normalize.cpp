#include "highmap/array.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/statistics.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Normalize, MinMax_RangeIs0To1)
{
  Array input = Array({{0.f, 5.f, 10.f}});
  Array out = normalized(input, NormalizationMethod::NM_MIN_MAX);

  EXPECT_NEAR(out.min(), 0.f, 1e-5);
  EXPECT_NEAR(out.max(), 1.f, 1e-5);
}

TEST(Normalize, MinMax_MonotonicityPreserved)
{
  Array input = Array({{1.f, 2.f, 3.f, 4.f}});
  Array out = normalized(input, NormalizationMethod::NM_MIN_MAX);

  for (int i = 1; i < input.shape.x; ++i)
    EXPECT_LE(out(i - 1, 0), out(i, 0));
}

TEST(Normalize, MinMax_ConstantField)
{
  Array input = Array({{5.f, 5.f, 5.f}});
  Array out = normalized(input, NormalizationMethod::NM_MIN_MAX);

  EXPECT_TRUE(assert_almost_equal(out, Array({{0.f, 0.f, 0.f}})));
}

TEST(Normalize, Standardize_MeanZero)
{
  Array input = Array({{1.f, 2.f, 3.f}});
  Array out = normalized(input, NormalizationMethod::NM_STANDARDIZE);

  float mean = out.mean();

  EXPECT_NEAR(mean, 0.f, 1e-5);
}

TEST(Normalize, Standardize_UnitVariance)
{
  Array input = Array({{1.f, 2.f, 3.f}});
  Array out = normalized(input, NormalizationMethod::NM_STANDARDIZE);

  float var = variance(out);

  EXPECT_NEAR(var, 1.f, 1e-4);
}

TEST(Normalize, Standardize_InvariantToShift)
{
  Array a = Array({{1.f, 2.f, 3.f}});
  Array b = Array({{101.f, 102.f, 103.f}});

  Array na = normalized(a, NormalizationMethod::NM_STANDARDIZE);
  Array nb = normalized(b, NormalizationMethod::NM_STANDARDIZE);

  EXPECT_TRUE(assert_almost_equal(na, nb));
}

TEST(Normalize, Robust_CenteringAroundMedian)
{
  Array input = Array({{1.f, 2.f, 100.f}});
  Array out = normalized(input, NormalizationMethod::NM_ROBUST);

  float median = out.median();

  EXPECT_NEAR(median, 0.f, 1e-5);
}

TEST(Normalize, Robust_ScaleInvariantToShift)
{
  Array a = Array({{1.f, 2.f, 3.f, 100.f}});
  Array b = Array({{10.f, 11.f, 12.f, 109.f}});

  Array na = normalized(a, NormalizationMethod::NM_ROBUST);
  Array nb = normalized(b, NormalizationMethod::NM_ROBUST);

  EXPECT_TRUE(assert_almost_equal(na, nb));
}

TEST(Normalize, Robust_HandlesSmallIQR)
{
  Array input = Array({{5.f, 5.f, 5.f, 6.f}});

  Array out = normalized(input, NormalizationMethod::NM_ROBUST);

  // should not explode
  for (int i = 0; i < input.size(); ++i)
    EXPECT_TRUE(std::isfinite(out.vector[i]));
}

TEST(Normalize, InvalidMethodThrows)
{
  Array input = Array({{1.f, 2.f, 3.f}});

  EXPECT_THROW(normalized(input, static_cast<NormalizationMethod>(999)),
               std::invalid_argument);
}
