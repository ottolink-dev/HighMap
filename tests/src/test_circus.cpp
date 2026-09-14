#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(Circus, BasicShapeAndBounds)
{
  glm::ivec2 shape = {128, 128};

  hmap::Array z = hmap::circus(shape);

  EXPECT_EQ(z.shape, shape);
  EXPECT_FALSE(std::isnan(z.min()));
  EXPECT_FALSE(std::isnan(z.max()));
  EXPECT_GE(z.min(), -1e-5f);
  EXPECT_LE(z.max(), 1.0f + 1e-5f);
}

TEST(Circus, CenterDepressionAndExitAndRidge)
{
  glm::ivec2 shape = {129, 129};
  glm::vec2  center = {0.5f, 0.5f};
  float      radius = 0.3f;
  float      exit_depth = 0.0f;
  float      center_depth = 0.1f;
  float      ridge_height = 1.0f;

  // angle = 0: exit point is towards +x (east), opposite ridge is towards -x
  // (west)
  hmap::Array z = hmap::circus(shape,
                               radius,
                               0.f, // angle
                               0.35f,
                               exit_depth,
                               center_depth,
                               ridge_height,
                               0.15f,
                               2.f,
                               nullptr,
                               center);

  // Center pixel (64, 64) should be close to center_depth
  EXPECT_NEAR(z(64, 64), center_depth, 0.02f);

  // Opposite ridge point around (center.x - radius, center.y) -> (64 - 0.3 *
  // 128, 64) ~= (26, 64)
  int ridge_i = 64 - static_cast<int>(radius * 128.f);
  EXPECT_NEAR(z(ridge_i, 64), ridge_height, 0.05f);

  // Exit point around (center.x + radius, center.y) -> (64 + 0.3 * 128, 64) ~=
  // (102, 64)
  int exit_i = 64 + static_cast<int>(radius * 128.f);
  EXPECT_NEAR(z(exit_i, 64), exit_depth, 0.05f);

  // Ridge should be significantly higher than depression and exit
  EXPECT_GT(z(ridge_i, 64), z(64, 64));
  EXPECT_GT(z(64, 64), z(exit_i, 64) - 1e-3f);
}

TEST(Circus, Rotation)
{
  glm::ivec2 shape = {129, 129};
  glm::vec2  center = {0.5f, 0.5f};
  float      radius = 0.3f;

  // Exit towards +x (angle = 0)
  hmap::Array z_east = hmap::circus(shape,
                                    radius,
                                    0.f,
                                    0.35f,
                                    0.f,
                                    0.1f,
                                    1.f,
                                    0.15f,
                                    2.f,
                                    nullptr,
                                    center);
  // Exit towards +y (angle = 90)
  hmap::Array z_north = hmap::circus(shape,
                                     radius,
                                     90.f,
                                     0.35f,
                                     0.f,
                                     0.1f,
                                     1.f,
                                     0.15f,
                                     2.f,
                                     nullptr,
                                     center);

  int offset = static_cast<int>(radius * 128.f);

  // East: Ridge at -x (64 - offset, 64), North: Ridge at -y (64, 64 - offset)
  EXPECT_NEAR(z_east(64 - offset, 64), z_north(64, 64 - offset), 1e-3f);

  // East: Exit at +x (64 + offset, 64), North: Exit at +y (64, 64 + offset)
  EXPECT_NEAR(z_east(64 + offset, 64), z_north(64, 64 + offset), 1e-3f);
}

TEST(Circus, MaskOutput)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array mask;

  hmap::Array z = hmap::circus(shape,
                               0.4f,
                               0.f,
                               0.35f,
                               0.f,
                               0.1f,
                               1.f,
                               0.15f,
                               2.f,
                               nullptr,
                               {0.5f, 0.5f},
                               {0.f, 1.f, 0.f, 1.f},
                               &mask);

  EXPECT_EQ(mask.shape, shape);
  EXPECT_NEAR(mask(32, 32), 1.0f, 1e-5f);
  EXPECT_GE(mask.min(), 0.f);
  EXPECT_LE(mask.max(), 1.f);
}

TEST(Circus, ContinuityAcrossRim)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  center = {0.5f, 0.5f};
  float      radius = 0.35f;

  hmap::Array z = hmap::circus(shape,
                               radius,
                               0.f,
                               0.35f,
                               0.0f,
                               0.1f,
                               1.0f,
                               0.15f,
                               2.0f,
                               nullptr,
                               center);

  // Across the entire map, adjacent pixel differences should be small
  // (continuous surface)
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 1; i < shape.x; ++i)
    {
      float diff_x = std::abs(z(i, j) - z(i - 1, j));
      EXPECT_LT(diff_x, 0.05f);
    }
  }

  for (int j = 1; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      float diff_y = std::abs(z(i, j) - z(i, j - 1));
      EXPECT_LT(diff_y, 0.05f);
    }
  }
}

TEST(Circus, NoiseModulation)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array noise = hmap::constant(shape, 0.2f);

  hmap::Array z_clean = hmap::circus(shape);
  hmap::Array z_noisy =
      hmap::circus(shape, 0.4f, 0.f, 0.35f, 0.f, 0.1f, 1.f, 0.15f, 2.f, &noise);

  EXPECT_FALSE(std::isnan(z_noisy.min()));
  EXPECT_NE(z_clean.mean(), z_noisy.mean());
}
