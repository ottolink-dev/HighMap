#include "highmap/dbg/assert.hpp"
#include "highmap/erosion.hpp"
#include "highmap/morphology.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Erosion, IdentityWhenRadiusZero)
{
  Array input = Array({{1, 2}, {3, 4}});
  Array cpu = erosion(input, 0);
  Array gpu = gpu::erosion(input, 0);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
  EXPECT_TRUE(assert_almost_equal(gpu, input));
}

TEST(Erosion, SingleMinimumPropagation)
{
  Array input = Array({{5, 5, 5}, {5, 1, 5}, {5, 5, 5}});

  // square kernel for GPU
  Array cpu = erosion(input, 1);
  Array expected_cpu = Array({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});

  // disk kernel for GPU
  Array gpu = gpu::erosion(input, 1);
  Array expected_gpu = Array({{5, 1, 5}, {1, 1, 1}, {5, 1, 5}});

  EXPECT_TRUE(assert_almost_equal(cpu, expected_cpu));
  EXPECT_TRUE(assert_almost_equal(gpu, expected_gpu));
}

TEST(Erosion, EdgeHandling)
{
  Array input = Array({{3, 3, 3}, {3, 1, 3}, {3, 3, 3}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  EXPECT_EQ(cpu(0, 0), 1);
  EXPECT_EQ(cpu(2, 2), 1);

  EXPECT_EQ(gpu(0, 0), 3);
  EXPECT_EQ(gpu(2, 2), 3);
}

TEST(Erosion, LargerRadius)
{
  Array input = Array({{9, 8, 7, 6}, {8, 5, 4, 7}, {7, 4, 3, 8}, {6, 7, 8, 9}});
  Array expected = Array(
      {{3, 3, 3, 3}, {3, 3, 3, 3}, {3, 3, 3, 3}, {3, 3, 3, 3}});

  Array cpu = erosion(input, 3);
  Array gpu = gpu::erosion(input, 3);

  EXPECT_TRUE(assert_almost_equal(cpu, expected));
  EXPECT_TRUE(assert_almost_equal(gpu, expected));
}

TEST(Erosion, MonotonicDecrease)
{
  Array input = Array({{5, 6}, {7, 8}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  for (int i = 0; i < input.shape.x; ++i)
    for (int j = 0; j < input.shape.y; ++j)
    {
      EXPECT_LE(cpu(i, j), input(i, j));
      EXPECT_LE(gpu(i, j), input(i, j));
    }
}

TEST(Erosion, FlatRegionUnchanged)
{
  Array input = Array({{2, 2, 2}, {2, 2, 2}, {2, 2, 2}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
  EXPECT_TRUE(assert_almost_equal(gpu, input));
}

TEST(Erosion, NonSquareArray)
{
  Array input = Array({{5, 4, 3, 2}});
  Array expected = Array({{4, 3, 2, 2}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, expected));
  EXPECT_TRUE(assert_almost_equal(gpu, expected));
}

TEST(ThermalConserve, MassPreserved)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  Array z = z0;
  gpu::thermal_conserve(z, 0.1f / shape.x, 100);

  EXPECT_NEAR(z.sum(), z0.sum(), 1e-2f);
}

TEST(ErosionFilter, BasicExecutionAndRidgeMap)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {2.f, 2.f};
  Array      z0 = noise_fbm(NoiseType::SIMPLEX2, shape, kw, 1234);

  Array z = z0;
  Array ridge_map(shape);

  gpu::erosion_filter(z,
                      /* scale */ 0.15f,
                      /* strength */ 0.22f,
                      /* gully_weight */ 0.5f,
                      /* detail */ 1.5f,
                      /* rounding */ {0.1f, 0.0f, 0.1f, 2.0f},
                      /* onset */ {1.25f, 1.25f, 2.8f, 1.5f},
                      /* assumed_slope */ {0.7f, 1.0f},
                      /* cell_scale */ 0.7f,
                      /* octaves */ 3,
                      /* gain */ 0.5f,
                      /* lacunarity */ 2.0f,
                      /* normalization */ 0.5f,
                      /* seed */ 1337,
                      /* p_fade_target */ nullptr,
                      /* p_ridge_map */ &ridge_map);

  // Check that heightmap was modified
  EXPECT_FALSE(assert_almost_equal(z, z0));

  // Check that ridge_map was populated
  EXPECT_GT(ridge_map.max(), ridge_map.min());

  // Determinism test with identical seed
  Array z_repeat = z0;
  gpu::erosion_filter(z_repeat,
                      0.15f,
                      0.22f,
                      0.5f,
                      1.5f,
                      {0.1f, 0.0f, 0.1f, 2.0f},
                      {1.25f, 1.25f, 2.8f, 1.5f},
                      {0.7f, 1.0f},
                      0.7f,
                      3,
                      0.5f,
                      2.0f,
                      0.5f,
                      1337);

  EXPECT_TRUE(assert_almost_equal(z, z_repeat));
}
