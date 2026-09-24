#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Canyonize, EmptyArray)
{
  Array z;
  canyonize(z);
  EXPECT_EQ(z.size(), 0);
}

TEST(Canyonize, PreservesShapeAndRange)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z, 0.f, 1.f);

  Array z_canyon = z;
  canyonize(z_canyon, 123, 5, 3.0f, 2.5f, 0.4f, 0.1f, 0.05f, 0.1f);

  EXPECT_EQ(z_canyon.shape.x, shape.x);
  EXPECT_EQ(z_canyon.shape.y, shape.y);
  EXPECT_GT(z_canyon.size(), 0);
  EXPECT_GE(z_canyon.min(), 0.09f); // bottom clamped smoothly around 0.1
  EXPECT_LE(z_canyon.max(), 1.05f);
}

TEST(Canyonize, FlattenBottom)
{
  glm::ivec2 shape = {32, 32};
  Array      z(shape, 0.05f); // Below clamp_min_val

  canyonize(z, 0, 4, 3.0f, 2.0f, 0.5f, 0.2f, 0.f);

  // Hard clamp at 0.2f should make all points 0.2f
  EXPECT_NEAR(z.min(), 0.2f, 1e-5);
  EXPECT_NEAR(z.max(), 0.2f, 1e-5);
}

TEST(Canyonize, MaskedCanyonize)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z, 0.f, 1.f);

  Array mask(shape, 0.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x / 2; ++i)
      mask(i, j) = 1.f;

  Array z_canyon = z;
  canyonize(z_canyon, 123, 6, &mask, 3.0f, 3.0f, 0.3f, 0.15f, 0.05f);

  // Right half should remain identical to original z
  for (int j = 0; j < shape.y; ++j)
    for (int i = shape.x / 2; i < shape.x; ++i)
      EXPECT_FLOAT_EQ(z_canyon(i, j), z(i, j));
}

TEST(Canyonize, LateralNoiseModulation)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z, 0.f, 1.f);

  Array noise = white(shape, 0.f, 0.05f, 99);

  Array z_with_noise = z;
  canyonize(z_with_noise, 42, 6, 3.0f, 3.0f, 0.3f, 0.1f, 0.05f, 0.1f, &noise);

  Array z_without_noise = z;
  canyonize(z_without_noise,
            42,
            6,
            3.0f,
            3.0f,
            0.3f,
            0.1f,
            0.05f,
            0.1f,
            nullptr);

  // They should differ because lateral noise modulates strata transitions
  bool diff_found = false;
  for (int i = 0; i < z.size(); ++i)
  {
    if (std::abs(z_with_noise.vector[i] - z_without_noise.vector[i]) > 1e-4f)
    {
      diff_found = true;
      break;
    }
  }
  EXPECT_TRUE(diff_found);
}
