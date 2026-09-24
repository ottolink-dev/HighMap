#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(QuadSurface, CornerValuesInterpolation)
{
  float c00 = 10.f;
  float c10 = 20.f;
  float c01 = 30.f;
  float c11 = 40.f;

  hmap::QuadSurfaceFunction f(c00, c10, c01, c11, {0.f, 1.f, 0.f, 1.f});

  // check corner values (xmin=0, ymin=0, xmax=1, ymax=1)
  EXPECT_NEAR(f.get_value(0.f, 0.f, 1.f), c00, 1e-5f);
  EXPECT_NEAR(f.get_value(1.f, 0.f, 1.f), c10, 1e-5f);
  EXPECT_NEAR(f.get_value(0.f, 1.f, 1.f), c01, 1e-5f);
  EXPECT_NEAR(f.get_value(1.f, 1.f, 1.f), c11, 1e-5f);

  // check center value (u = 0.5, v = 0.5)
  EXPECT_NEAR(f.get_value(0.5f, 0.5f, 1.f),
              0.25f * (c00 + c10 + c01 + c11),
              1e-5f);
}

TEST(QuadSurface, ArrayEvaluation)
{
  glm::ivec2 shape = {64, 64};

  float c00 = 0.f;
  float c10 = 1.f;
  float c01 = 1.f;
  float c11 = 0.f;

  hmap::Array z = hmap::quad_surface(shape, c00, c10, c01, c11);

  EXPECT_EQ(z.shape, shape);
  EXPECT_GE(z.min(), 0.f - 1e-5f);
  EXPECT_LE(z.max(), 1.f + 1e-5f);
}

TEST(QuadSurface, FlatSurface)
{
  glm::ivec2 shape = {32, 32};
  float      val = 5.f;

  hmap::Array z = hmap::quad_surface(shape, val, val, val, val);

  EXPECT_NEAR(z.min(), val, 1e-6f);
  EXPECT_NEAR(z.max(), val, 1e-6f);
}

TEST(QuadSurface, FunctionDirectEvaluation)
{
  hmap::QuadSurfaceFunction f(0.f, 10.f, 20.f, 30.f, {0.f, 2.f, 0.f, 4.f});

  // Evaluate at corners
  EXPECT_NEAR(f.get_value(0.f, 0.f, 1.f), 0.f, 1e-6f);
  EXPECT_NEAR(f.get_value(2.f, 0.f, 1.f), 10.f, 1e-6f);
  EXPECT_NEAR(f.get_value(0.f, 4.f, 1.f), 20.f, 1e-6f);
  EXPECT_NEAR(f.get_value(2.f, 4.f, 1.f), 30.f, 1e-6f);

  // Evaluate at midpoint (x=1.0, y=2.0)
  EXPECT_NEAR(f.get_value(1.f, 2.f, 1.f), 15.f, 1e-6f);

  // With ctrl_param multiplier
  EXPECT_NEAR(f.get_value(1.f, 2.f, 2.f), 30.f, 1e-6f);
}
