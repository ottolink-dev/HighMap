#include "highmap/curvature.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;
TEST(LevelSetCurvature, CpuGpuFlatField)
{
  Array z{
      {1.f, 1.f, 1.f},
      {1.f, 1.f, 1.f},
      {1.f, 1.f, 1.f},
  };

  int ir = 0;

  Array c0 = level_set_curvature(z, ir);
  Array c1 = gpu::level_set_curvature(z, ir);

  bool ret = assert_almost_equal(c0, c1, 1e-5f);
  EXPECT_TRUE(ret);

  Array expected(z.shape, 0.f);

  ret = assert_almost_equal(c0, expected, 1e-5f);
  EXPECT_TRUE(ret);
}

TEST(LevelSetCurvature, CpuGpuLinearRampX)
{
  Array z{
      {0.f, 1.f, 2.f, 3.f},
      {0.f, 1.f, 2.f, 3.f},
      {0.f, 1.f, 2.f, 3.f},
  };

  int ir = 0;

  Array c0 = level_set_curvature(z, ir);
  Array c1 = gpu::level_set_curvature(z, ir);

  bool ret = assert_almost_equal(c0, c1, 1e-4f);
  EXPECT_TRUE(ret);

  Array expected(z.shape, 0.f);

  ret = assert_almost_equal(c0, expected, 1e-4f);
  EXPECT_TRUE(ret);
}

TEST(LevelSetCurvature, CpuGpuParaboloid)
{
  Array z{
      {2.f, 1.f, 2.f},
      {1.f, 0.f, 1.f},
      {2.f, 1.f, 2.f},
  };

  int ir = 0;

  Array c0 = level_set_curvature(z, ir);
  Array c1 = gpu::level_set_curvature(z, ir);

  bool ret = assert_almost_equal(c0, c1, 1e-4f);
  EXPECT_TRUE(ret);
}

TEST(LevelSetCurvature, CpuGpuWithPrefilter)
{
  Array z{
      {0.f, 1.f, 0.f, 1.f},
      {1.f, 0.f, 1.f, 0.f},
      {0.f, 1.f, 0.f, 1.f},
      {1.f, 0.f, 1.f, 0.f},
  };

  int ir = 2;

  Array c0 = level_set_curvature(z, ir);
  Array c1 = gpu::level_set_curvature(z, ir);

  bool ret = assert_almost_equal(c0, c1, 1e-3f);
  EXPECT_TRUE(ret);
}

TEST(LevelSetCurvature, CpuGpuRandomField)
{
  Array z = white({64, 64}, -1.f, 1.f, 42);

  int ir = 1;

  Array c0 = level_set_curvature(z, ir);
  Array c1 = gpu::level_set_curvature(z, ir);

  bool ret = assert_almost_equal(c0, c1, 1e-3f);
  EXPECT_TRUE(ret);
}

TEST(LevelSetCurvature, OutputShape)
{
  Array z = white({32, 16}, 0.f, 1.f, 123);

  int ir = 0;

  Array c0 = level_set_curvature(z, ir);

  EXPECT_EQ(c0.shape.x, z.shape.x);
  EXPECT_EQ(c0.shape.y, z.shape.y);
}
