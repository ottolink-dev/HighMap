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

TEST(ThermalGPU, VariantsRunAndModify)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  float      talus = 0.5f / shape.x;
  Array      talus_map(shape, talus);

  // thermal
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal(z, talus_map, 5, nullptr, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_auto_bedrock
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal_auto_bedrock(z, talus_map, 5, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_flatten
  {
    Array z = z0;
    gpu::thermal_flatten(z, talus_map, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_inflate
  {
    Array z = z0;
    gpu::thermal_inflate(z, talus_map, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_olsen
  {
    Array z = z0;
    gpu::thermal_olsen(z, talus_map, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_rib
  {
    Array z = z0;
    gpu::thermal_rib(z, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_ridge
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal_ridge(z, talus_map, 5, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_schott
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal_schott(z, talus_map, 5, 0.2f, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_scree
  {
    Array z = z0;
    Array zmax(shape, 1.f);
    Array dep(shape);
    gpu::thermal_scree(z, talus_map, zmax, 5, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }
}
