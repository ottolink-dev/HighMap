#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(AdaptiveRelief, ConstantArrayInvariant)
{
  const glm::ivec2 shape = {64, 64};
  Array            z = Array(shape, 5.0f);
  Array            z_orig = z;

  adaptive_relief(z, 1.0f, 1.0f);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      EXPECT_NEAR(z(i, j), z_orig(i, j), 1e-5f);
}

TEST(AdaptiveRelief, BoundedRangeClamp)
{
  const glm::ivec2 shape = {128, 128};
  const glm::vec2  kw = {4.f, 4.f};
  const int        seed = 123;

  Array z = noise_fbm(NoiseType::PERLIN, shape, kw, seed);
  float min_orig = z.min();
  float max_orig = z.max();

  adaptive_relief(z, 2.0f, 1.0f); // strict clamp

  EXPECT_GE(z.min(), min_orig - 1e-5f);
  EXPECT_LE(z.max(), max_orig + 1e-5f);
}

TEST(AdaptiveRelief, MaskSupport)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const int        seed = 10;

  Array z = noise_fbm(NoiseType::PERLIN, shape, kw, seed);
  Array z_orig = z;

  Array mask = Array(shape, 0.0f); // 0 everywhere -> no change expected
  adaptive_relief(z, &mask, 1.5f, 1.0f);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      EXPECT_NEAR(z(i, j), z_orig(i, j), 1e-5f);
}

TEST(AdaptiveRelief, CpuGpuConsistency)
{
  const glm::ivec2 shape = {256, 256};
  const glm::vec2  kw = {4.f, 4.f};
  const int        seed = 42;

  Array z = noise_fbm(NoiseType::PERLIN, shape, kw, seed);
  remap(z, 0.f, 1.f);

  Array z_cpu = z;
  Array z_gpu = z;

  adaptive_relief(z_cpu, 1.5f, 1.0f);
  gpu::adaptive_relief(z_gpu, 1.5f, 1.0f);

  bool ret = assert_almost_equal(z_cpu, z_gpu, 1e-4f);
  EXPECT_TRUE(ret);
}
