#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(RecastCliff, CliffMaskFlatTerrainIsZero)
{
  Array z = Array(glm::ivec2(64, 64), 0.5f);
  Array cliff_mask;

  gpu::recast_cliff(z, 0.1f, 4, 0.2f, 2.f, 50, &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, 64);
  EXPECT_EQ(cliff_mask.shape.y, 64);
  EXPECT_NEAR(cliff_mask.max(), 0.f, 1e-6);
  EXPECT_NEAR(cliff_mask.min(), 0.f, 1e-6);
}

TEST(RecastCliff, CliffMaskSteepTerrainIsPositive)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float gain = 2.f;

  gpu::recast_cliff(z_copy, talus, ir, amplitude, gain, 50, &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);
  EXPECT_GT(cliff_mask.max(), 0.f);
}

TEST(RecastCliff, CliffMaskWithInputMask)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array mask = Array(shape, 0.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x / 2; ++i)
      mask(i, j) = 1.f;

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float gain = 2.f;

  gpu::recast_cliff(z_copy, talus, ir, amplitude, &mask, gain, 50, &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);

  // right half should be zero due to mask
  for (int j = 0; j < shape.y; ++j)
    for (int i = shape.x / 2; i < shape.x; ++i)
      EXPECT_NEAR(cliff_mask(i, j), 0.f, 1e-6);
}

TEST(RecastCliffDirectional, CliffMaskPositive)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float angle = 45.f;
  float gain = 2.f;

  gpu::recast_cliff_directional(z_copy,
                                talus,
                                ir,
                                amplitude,
                                angle,
                                gain,
                                50,
                                &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);
  EXPECT_GT(cliff_mask.max(), 0.f);
}

TEST(RecastCliffDirectional, CliffMaskWithInputMask)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array mask = Array(shape, 0.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x / 2; ++i)
      mask(i, j) = 1.f;

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float angle = 45.f;
  float gain = 2.f;

  gpu::recast_cliff_directional(z_copy,
                                talus,
                                ir,
                                amplitude,
                                angle,
                                &mask,
                                gain,
                                50,
                                &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);

  // right half should be zero due to mask
  for (int j = 0; j < shape.y; ++j)
    for (int i = shape.x / 2; i < shape.x; ++i)
      EXPECT_NEAR(cliff_mask(i, j), 0.f, 1e-6);
}

TEST(RecastCliffDirectional, CliffMaskWithVariableAngle)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array angle_field = noise_fbm(NoiseType::PERLIN,
                                shape,
                                glm::vec2(1.f, 1.f),
                                10);
  remap(angle_field, 0.f, 360.f);

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float gain = 2.f;

  gpu::recast_cliff_directional(z_copy,
                                talus,
                                ir,
                                amplitude,
                                angle_field,
                                gain,
                                50,
                                &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);
  EXPECT_GT(cliff_mask.max(), 0.f);
}
