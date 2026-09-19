#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(ValleyHead, HeadAndUpstreamZero)
{
  glm::ivec2 shape = {128, 128};
  glm::vec2  center = {0.5f, 0.5f};

  // angle = 0: valley flows toward +y (from y = 0.5 to y = 1.0)
  hmap::Array z = hmap::valley_head(shape,
                                    0.f,   // angle
                                    0.2f,  // max_depth
                                    0.05f, // min_width
                                    0.2f,  // max_width
                                    0.4f,  // depth_length
                                    0.4f,  // width_length
                                    2.f,   // development_power
                                    1.3f,  // profile_power
                                    nullptr,
                                    nullptr,
                                    center);

  // at valley head (x = 0.5, y = 0.5 -> grid index (64, 64))
  EXPECT_NEAR(z(64, 64), 0.f, 1e-5f);

  // upstream (y < 0.5, i.e., rows j < 64): must be strictly zero everywhere
  for (int j = 0; j <= 64; ++j)
    for (int i = 0; i < shape.x; ++i)
      EXPECT_NEAR(z(i, j), 0.f, 1e-6f);
}

TEST(ValleyHead, DownstreamDepthIncreases)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  center = {0.5f, 0.2f};
  float      max_depth = 0.4f;

  hmap::Array z = hmap::valley_head(shape,
                                    0.f,
                                    max_depth,
                                    0.05f,
                                    0.2f,
                                    0.3f,
                                    0.3f,
                                    2.f,
                                    1.3f,
                                    nullptr,
                                    nullptr,
                                    center);

  // sample along centerline x = 0.5 (column 128) moving downstream (y
  // increasing)
  int   center_col = 128;
  float prev_incision = 0.f;

  for (int j = 60; j < 240; j += 20)
  {
    float incision = -z(center_col, j); // incision depth (positive)
    EXPECT_GE(incision, prev_incision - 1e-5f);
    EXPECT_LE(incision, max_depth + 1e-5f);
    prev_incision = incision;
  }
}

TEST(ValleyHead, WidthIncreasesDownstream)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  center = {0.5f, 0.2f};

  hmap::Array z = hmap::valley_head(shape,
                                    0.f,
                                    0.3f,
                                    0.04f,
                                    0.25f,
                                    0.3f,
                                    0.3f,
                                    2.f,
                                    1.3f,
                                    nullptr,
                                    nullptr,
                                    center);

  // measure half-width (extent where z < 0) at near-head row vs downstream row
  auto get_half_width_pixels = [&](int row) -> int
  {
    int center_col = 128;
    int half_width = 0;
    for (int i = center_col; i < shape.x; ++i)
    {
      if (z(i, row) < -1e-5f) half_width = i - center_col;
    }
    return half_width;
  };

  int width_near_head = get_half_width_pixels(70);
  int width_downstream = get_half_width_pixels(220);

  EXPECT_GT(width_downstream, width_near_head);
}

TEST(ValleyHead, SymmetricCrossSection)
{
  glm::ivec2 shape = {128, 128};
  glm::vec2  center = {0.5f, 0.2f};

  hmap::Array z = hmap::valley_head(shape,
                                    0.f,
                                    0.3f,
                                    0.05f,
                                    0.2f,
                                    0.3f,
                                    0.3f,
                                    2.f,
                                    1.5f,
                                    nullptr,
                                    nullptr,
                                    center);

  int center_col = 64;
  for (int j = 60; j < 120; j += 10)
  {
    for (int offset = 1; offset < 30; ++offset)
    {
      EXPECT_NEAR(z(center_col - offset, j), z(center_col + offset, j), 1e-5f);
    }
  }
}

TEST(ValleyHead, OutsideInfluenceZero)
{
  glm::ivec2 shape = {128, 128};
  glm::vec2  center = {0.5f, 0.2f};

  hmap::Array z = hmap::valley_head(shape,
                                    0.f,
                                    0.3f,
                                    0.02f,
                                    0.1f,
                                    0.3f,
                                    0.3f,
                                    2.f,
                                    1.3f,
                                    nullptr,
                                    nullptr,
                                    center);

  // far left and far right columns should be zero everywhere
  for (int j = 0; j < shape.y; ++j)
  {
    EXPECT_NEAR(z(5, j), 0.f, 1e-6f);
    EXPECT_NEAR(z(120, j), 0.f, 1e-6f);
  }
}

TEST(ValleyHead, LateralContinuity)
{
  glm::ivec2 shape = {512, 512};
  glm::vec2  center = {0.5f, 0.2f};

  hmap::Array z = hmap::valley_head(shape,
                                    0.f,
                                    0.5f,
                                    0.05f,
                                    0.2f,
                                    0.4f,
                                    0.4f,
                                    2.f,
                                    1.3f,
                                    nullptr,
                                    nullptr,
                                    center);

  // across cross-sections, differences between adjacent pixels should be small
  int row = 350;
  for (int i = 1; i < shape.x; ++i)
  {
    float diff = std::abs(z(i, row) - z(i - 1, row));
    EXPECT_LT(diff, 0.02f);
  }
}

TEST(ValleyHead, ParameterExtremes)
{
  glm::ivec2 shape = {64, 64};

  // min_width == max_width
  hmap::Array z1 = hmap::valley_head(shape, 0.f, 0.2f, 0.1f, 0.1f);
  EXPECT_LE(z1.min(), 0.f);
  EXPECT_GE(z1.max(), 0.f);

  // profile_power == 1
  hmap::Array z2 =
      hmap::valley_head(shape, 0.f, 0.2f, 0.05f, 0.15f, 0.5f, 0.5f, 2.f, 1.f);
  EXPECT_LE(z2.min(), 0.f);

  // large profile_power
  hmap::Array z3 =
      hmap::valley_head(shape, 0.f, 0.2f, 0.05f, 0.15f, 0.5f, 0.5f, 2.f, 5.f);
  EXPECT_LE(z3.min(), 0.f);

  // small length scales
  hmap::Array z4 =
      hmap::valley_head(shape, 0.f, 0.2f, 0.05f, 0.15f, 0.01f, 0.01f);
  EXPECT_LE(z4.min(), 0.f);
}

TEST(ValleyHead, OrientationRotation)
{
  glm::ivec2 shape = {128, 128};
  glm::vec2  center = {0.5f, 0.5f};

  // 0 degrees (+y) vs 90 degrees (+x)
  hmap::Array z_north = hmap::valley_head(shape,
                                          0.f,
                                          0.3f,
                                          0.05f,
                                          0.2f,
                                          0.3f,
                                          0.3f,
                                          2.f,
                                          1.3f,
                                          nullptr,
                                          nullptr,
                                          center);
  hmap::Array z_east = hmap::valley_head(shape,
                                         90.f,
                                         0.3f,
                                         0.05f,
                                         0.2f,
                                         0.3f,
                                         0.3f,
                                         2.f,
                                         1.3f,
                                         nullptr,
                                         nullptr,
                                         center);

  // down valley by delta = +20 pixels from center (64, 64)
  // for north: (64, 84)
  // for east:  (84, 64)
  EXPECT_NEAR(z_north(64, 84), z_east(84, 64), 1e-4f);

  // cross valley at delta y = +20, delta x = +10
  // for north: (74, 84)
  // for east: (84, 54) (since x_local = -dy for 90 deg)
  EXPECT_NEAR(z_north(74, 84), z_east(84, 54), 1e-4f);
}

TEST(ValleyHead, CrossSectionNoise)
{
  glm::ivec2 shape = {128, 128};
  glm::vec2  center = {0.5f, 0.2f};

  hmap::Array noise = hmap::constant(
      shape,
      0.5f); // positive offset shifts the valley centerline

  hmap::Array z_clean = hmap::valley_head(shape,
                                          0.f,
                                          0.3f,
                                          0.05f,
                                          0.2f,
                                          0.3f,
                                          0.3f,
                                          2.f,
                                          1.3f,
                                          nullptr,
                                          nullptr,
                                          center);
  hmap::Array z_noisy = hmap::valley_head(shape,
                                          0.f,
                                          0.3f,
                                          0.05f,
                                          0.2f,
                                          0.3f,
                                          0.3f,
                                          2.f,
                                          1.3f,
                                          &noise,
                                          nullptr,
                                          center);

  // at centerline (64, 80), z_noisy should differ from z_clean because center
  // is shifted
  EXPECT_NE(z_clean(64, 80), z_noisy(64, 80));
}
