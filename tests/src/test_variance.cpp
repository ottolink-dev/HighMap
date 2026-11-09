#include "highmap/dbg/assert.hpp"
#include "highmap/statistics.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Variance, ConstantArray)
{
  Array input = Array({{2, 2, 2}, {2, 2, 2}});

  float var = variance(input, nullptr);

  EXPECT_FLOAT_EQ(var, 0.f);
}

TEST(Variance, KnownValues)
{
  // mean = 2.5
  // variance = ((1.5^2 + 0.5^2 + 0.5^2 + 1.5^2) / 4) = 1.25
  Array input = Array({{1, 2}, {3, 4}});

  float var = variance(input, nullptr);

  EXPECT_NEAR(var, 1.25f, 1e-6f);
}

TEST(Variance, WithProvidedMean)
{
  Array input = Array({{1, 2}, {3, 4}});

  float mean = input.mean();

  float var1 = variance(input, nullptr);
  float var2 = variance(input, &mean);

  EXPECT_NEAR(var1, var2, 1e-6f);
}

TEST(Variance, WrongMeanChangesResult)
{
  Array input = Array({{1, 2}, {3, 4}});

  float wrong_mean = 0.f;

  float var_correct = variance(input, nullptr);
  float var_wrong = variance(input, &wrong_mean);

  EXPECT_NE(var_correct, var_wrong);
}

TEST(Variance, InvariantToOffset)
{
  Array input = Array({{1, 2}, {3, 4}});

  Array shifted = input + 10.f;

  float v1 = variance(input, nullptr);
  float v2 = variance(shifted, nullptr);

  EXPECT_NEAR(v1, v2, 1e-6f);
}

TEST(Variance, ScalingProperty)
{
  Array input = Array({{1, 2}, {3, 4}});

  Array scaled = input * 3.f;

  float v1 = variance(input, nullptr);
  float v2 = variance(scaled, nullptr);

  EXPECT_NEAR(v2, 9.f * v1, 1e-6f); // variance scales with square
}

TEST(Variance, NonSquareArray)
{
  Array input = Array({{1, 2, 3, 4}});

  float var = variance(input, nullptr);

  EXPECT_GT(var, 0.f);
}

TEST(ArrayStd, MeanOutputPointer)
{
  Array input = Array({{1.f, 2.f}, {3.f, 4.f}});

  float mean_out = 0.f;
  float std_val = input.std(&mean_out);

  EXPECT_NEAR(mean_out, 2.5f, 1e-6f);
  EXPECT_NEAR(std_val, std::sqrt(1.25f), 1e-6f);
}
