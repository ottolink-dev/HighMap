#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(SnowSimulation, ZeroIterationsReturnsEmpty)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array fall_map(shape, 1.f);
  Array melting_map(shape, 0.f);
  Array talus(shape, 2.5f / shape.x);

  Array s =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 0, 0.5f);

  EXPECT_EQ(s.shape, shape);
  EXPECT_FLOAT_EQ(s.max(), 0.f);
  EXPECT_FLOAT_EQ(s.min(), 0.f);
}

TEST(SnowSimulation, IterationParityAndDeterminism)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array fall_map(shape, 1.f);
  Array melting_map(shape, 0.f);
  Array talus(shape, 2.5f / shape.x);

  Array s_even =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 50, 0.5f);
  Array s_odd =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 51, 0.5f);
  Array s_even_again =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 50, 0.5f);

  // deterministic within float precision
  EXPECT_NEAR(s_even.sum(), s_even_again.sum(), 1e-2f);

  // evolution check
  const float sum_even = s_even.sum();
  const float sum_odd = s_odd.sum();
  EXPECT_GT(sum_even, 0.f);
  EXPECT_NEAR(sum_odd, sum_even, 0.05f * sum_even);
  EXPECT_GE(s_odd.min(), 0.f);
}

TEST(SnowSimulation, SnowAccumulationWithoutMelting)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 7);
  remap(z);

  Array fall_map(shape, 1.f);
  Array melting_map(shape, 0.f);
  Array talus(shape, 2.5f / shape.x);

  Array s =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 100, 0.5f);

  EXPECT_GT(s.sum(), 0.f);
  EXPECT_GE(s.min(), 0.f);
}

TEST(SnowSimulation, MeltingReducesSnowVolume)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array fall_map(shape, 1.f);
  Array no_melt(shape, 0.f);
  Array melt(shape, 1.f);
  Array talus(shape, 2.5f / shape.x);

  Array s_no_melt =
      gpu::snow_simulation(z, 0.1f, fall_map, no_melt, talus, 100, 0.5f);
  Array s_melt =
      gpu::snow_simulation(z, 0.1f, fall_map, melt, talus, 100, 0.5f);

  EXPECT_GT(s_no_melt.sum(), s_melt.sum());
  EXPECT_GE(s_melt.min(), 0.f);
}
