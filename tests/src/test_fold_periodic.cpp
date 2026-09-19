#include "highmap.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(FoldPeriodic, EmptyArray)
{
  Array array;
  fold_periodic(array, PhasorProfile::PP_COSINE_STD);
  EXPECT_TRUE(array.vector.empty());
}

TEST(FoldPeriodic, SingleIterationCosineStd)
{
  Array array = Array({{0.f, 0.5f, 1.f}});
  // With vmin=0, vmax=1, iterations=1, it=0 -> freq=1:
  // t = 0   -> phi = -pi -> cos(-pi) = -1 -> -cos(-pi) = 1
  // t = 0.5 -> phi = 0   -> cos(0) = 1   -> -cos(0) = -1
  // t = 1   -> phi = pi  -> cos(pi) = -1 -> -cos(pi) = 1
  fold_periodic(array, PhasorProfile::PP_COSINE_STD, 0.f, 1.f, 1);

  EXPECT_NEAR(array(0, 0), 1.f, 1e-5f);
  EXPECT_NEAR(array(1, 0), -1.f, 1e-5f);
  EXPECT_NEAR(array(2, 0), 1.f, 1e-5f);
}

TEST(FoldPeriodic, SingleIterationTriangle)
{
  Array array = Array({{0.f, 0.5f, 1.f}});
  // PP_TRIANGLE:
  // t = 0   -> phi = -pi -> t_tri = -1 -> 4 * abs(-0.5 - floor(0)) - 1 = 4 *
  // 0.5 - 1 = 1 t = 0.5 -> phi = 0   -> t_tri = 0  -> 4 * abs(0 - floor(0.5)) -
  // 1 = 4 * 0 - 1 = -1 t = 1   -> phi = pi  -> t_tri = 1  -> 4 * abs(0.5 -
  // floor(1)) - 1 = 4 * 0.5 - 1 = 1
  fold_periodic(array, PhasorProfile::PP_TRIANGLE, 0.f, 1.f, 1);

  EXPECT_NEAR(array(0, 0), 1.f, 1e-5f);
  EXPECT_NEAR(array(1, 0), -1.f, 1e-5f);
  EXPECT_NEAR(array(2, 0), 1.f, 1e-5f);
}

TEST(FoldPeriodic, DefaultRangeOverload)
{
  Array array = Array({{10.f, 15.f, 20.f}});
  fold_periodic(array, PhasorProfile::PP_COSINE_STD, 1);

  EXPECT_NEAR(array(0, 0), 1.f, 1e-5f);
  EXPECT_NEAR(array(1, 0), -1.f, 1e-5f);
  EXPECT_NEAR(array(2, 0), 1.f, 1e-5f);
}

TEST(FoldPeriodic, DegenerateConstantRange)
{
  Array array = Array({{5.f, 5.f, 5.f}});
  fold_periodic(array, PhasorProfile::PP_COSINE_STD, 5.f, 5.f, 1);

  EXPECT_NEAR(array(0, 0), 0.f, 1e-5f);
  EXPECT_NEAR(array(1, 0), 0.f, 1e-5f);
  EXPECT_NEAR(array(2, 0), 0.f, 1e-5f);
}

TEST(FoldPeriodic, MultipleIterations)
{
  glm::ivec2 shape = {64, 64};
  Array      z = white(shape, 0.f, 1.f, 42);

  Array z1 = z;
  fold_periodic(z1, PhasorProfile::PP_COSINE_STD, 0.f, 1.f, 3);

  EXPECT_GE(z1.min(), -1.01f);
  EXPECT_LE(z1.max(), 1.01f);
}
