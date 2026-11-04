#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(SnowMantle, EmptyInputReturnsEmpty)
{
  Array z;
  Array s = gpu::snow_mantle(z);
  EXPECT_EQ(s.size(), 0);
}

TEST(SnowMantle, BasicAccumulation)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z, 0.f, 1.f);

  Array s = gpu::snow_mantle(z,
                             0.1f,                  // strength
                             0.5f,                  // snowline
                             0.1f,                  // falloff
                             8,                     // prefilter_ir
                             0.8f,                  // talus_coverage
                             0.6f / float(shape.x), // talus_threshold
                             0.0f,                  // directional
                             0.0f,                  // wind_direction
                             0.0f,                  // peaks
                             0.0f,                  // melt
                             180.f,                 // sun_azimuth
                             45.f,                  // sun_zenith
                             0,                     // thermal_iterations
                             0.15f,                 // thermal_talus
                             0.5f);                 // coverage

  EXPECT_EQ(s.shape, shape);
  EXPECT_GE(s.min(), 0.f);
  EXPECT_GT(s.max(), 0.f);

  // snow should be concentrated above snowline - falloff = 0.4
  float low_alt_snow = 0.f;
  float high_alt_snow = 0.f;

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      if (z(i, j) < 0.3f)
        low_alt_snow += s(i, j);
      else if (z(i, j) > 0.7f)
        high_alt_snow += s(i, j);
    }
  }

  EXPECT_LT(low_alt_snow, high_alt_snow);
}

TEST(SnowMantle, FlowPreservesNonNegativityAndSmooths)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z, 0.f, 1.f);

  Array s_noflow = gpu::snow_mantle(z,
                                    0.1f,
                                    0.4f,
                                    0.1f,
                                    8,
                                    0.5f,
                                    0.5f / float(shape.x),
                                    0.f,
                                    0.f,
                                    0.f,
                                    0.f,
                                    180.f,
                                    45.f,
                                    0);
  Array s_flow = gpu::snow_mantle(z,
                                  0.1f,
                                  0.4f,
                                  0.1f,
                                  8,
                                  0.5f,
                                  0.5f / float(shape.x),
                                  0.f,
                                  0.f,
                                  0.f,
                                  0.f,
                                  180.f,
                                  45.f,
                                  20);

  EXPECT_GE(s_noflow.min(), 0.f);
  EXPECT_GE(s_flow.min(), 0.f);
  EXPECT_GT(s_flow.sum(), 0.f);
}

TEST(SnowMantle, MeltingReducesSnowVolume)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z, 0.f, 1.f);

  Array s_nomelt = gpu::snow_mantle(z,
                                    0.1f,
                                    0.3f,
                                    0.1f,
                                    8,
                                    0.5f,
                                    0.5f / float(shape.x),
                                    0.f,
                                    0.f,
                                    0.f,
                                    0.f,
                                    180.f,
                                    45.f,
                                    0);
  Array s_melt = gpu::snow_mantle(z,
                                  0.1f,
                                  0.3f,
                                  0.1f,
                                  8,
                                  0.5f,
                                  0.5f / float(shape.x),
                                  0.f,
                                  0.f,
                                  0.f,
                                  0.8f,
                                  180.f,
                                  45.f,
                                  0);

  EXPECT_GT(s_nomelt.sum(), s_melt.sum());
}

TEST(SnowMantle, DirectionalBias)
{
  gpu::init_opencl();

  // Create an East-West ramp: slopes face West or East
  const glm::ivec2 shape = {64, 64};
  Array            z(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      z(i, j) = float(i) / float(shape.x);

  // Wind blowing towards East (0 deg) vs West (180 deg)
  Array s_east = gpu::snow_mantle(z,
                                  0.1f,
                                  0.2f,
                                  0.1f,
                                  8,
                                  0.f,
                                  1.f,
                                  0.5f,
                                  0.f,
                                  0.f,
                                  0.f,
                                  180.f,
                                  45.f,
                                  0);
  Array s_west = gpu::snow_mantle(z,
                                  0.1f,
                                  0.2f,
                                  0.1f,
                                  8,
                                  0.f,
                                  1.f,
                                  0.5f,
                                  180.f,
                                  0.f,
                                  0.f,
                                  180.f,
                                  45.f,
                                  0);

  // Normal dot wind directions differ
  EXPECT_NE(s_east.sum(), s_west.sum());
}

TEST(SnowMantle, NoiseModulation)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z, 0.f, 1.f);

  Array noise(shape, 0.5f);

  Array s_plain = gpu::snow_mantle(z,
                                   0.1f,
                                   0.3f,
                                   0.1f,
                                   8,
                                   0.5f,
                                   0.5f / float(shape.x),
                                   0.f,
                                   0.f,
                                   0.f,
                                   0.f,
                                   180.f,
                                   45.f,
                                   0,
                                   0.15f,
                                   0.f);
  Array s_noisy = gpu::snow_mantle(z,
                                   0.1f,
                                   0.3f,
                                   0.1f,
                                   8,
                                   0.5f,
                                   0.5f / float(shape.x),
                                   0.f,
                                   0.f,
                                   0.f,
                                   0.f,
                                   180.f,
                                   45.f,
                                   0,
                                   0.15f,
                                   0.f,
                                   nullptr,
                                   nullptr,
                                   &noise);

  EXPECT_GT(s_plain.sum(), 0.f);
  EXPECT_GT(s_noisy.sum(), 0.f);
  EXPECT_NEAR(s_noisy.sum(), 0.5f * s_plain.sum(), 1e-4f);
}
