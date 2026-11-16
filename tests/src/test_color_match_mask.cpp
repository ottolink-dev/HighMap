#include "highmap/colorize.hpp"
#include "highmap/internal/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(ColorMatchMask, ExactMatch)
{
  Array r = Array({{0.2f}});
  Array g = Array({{0.4f}});
  Array b = Array({{0.6f}});

  glm::vec3 target(0.2f, 0.4f, 0.6f);

  Array mask = color_match_mask(r, g, b, target, 0.f);

  EXPECT_EQ(mask(0, 0), 1.f);
}

TEST(ColorMatchMask, WithinToleranceSmallDeviation)
{
  Array r = Array({{0.2f}});
  Array g = Array({{0.4f}});
  Array b = Array({{0.6f}});

  glm::vec3 target(0.21f, 0.41f, 0.61f);

  Array mask = color_match_mask(r, g, b, target, 0.05f);

  EXPECT_EQ(mask(0, 0), 1.f);
}

TEST(ColorMatchMask, OutsideTolerance)
{
  Array r = Array({{0.2f}});
  Array g = Array({{0.4f}});
  Array b = Array({{0.6f}});

  glm::vec3 target(0.5f, 0.5f, 0.5f);

  Array mask = color_match_mask(r, g, b, target, 0.1f);

  EXPECT_EQ(mask(0, 0), 0.f);
}

TEST(ColorMatchMask, ExactlyOnToleranceBoundary)
{
  // distance = 0.1
  Array r = Array({{0.0f}});
  Array g = Array({{0.0f}});
  Array b = Array({{0.0f}});

  glm::vec3 target(0.1f, 0.f, 0.f);

  Array mask = color_match_mask(r, g, b, target, 0.1f);

  EXPECT_EQ(mask(0, 0), 1.f);
}

TEST(ColorMatchMask, MultiPixel)
{
  Array r = Array({{0.2f, 0.9f}, {0.2f, 0.9f}});
  Array g = Array({{0.4f, 0.9f}, {0.4f, 0.9f}});
  Array b = Array({{0.6f, 0.9f}, {0.6f, 0.9f}});

  glm::vec3 target(0.2f, 0.4f, 0.6f);

  Array mask = color_match_mask(r, g, b, target, 0.05f);

  EXPECT_EQ(mask(0, 0), 1.f);
  EXPECT_EQ(mask(0, 1), 0.f);
  EXPECT_EQ(mask(1, 0), 1.f);
  EXPECT_EQ(mask(1, 1), 0.f);
}

TEST(ColorMatchMask, VerySmallTolerance)
{
  Array r = Array({{0.2f}});
  Array g = Array({{0.4f}});
  Array b = Array({{0.6f}});

  glm::vec3 target(0.2001f, 0.4f, 0.6f);

  Array mask = color_match_mask(r, g, b, target, 1e-5f);

  EXPECT_EQ(mask(0, 0), 0.f);
}

TEST(ColorMatchMask, MaxRangeDifference)
{
  Array r = Array({{0.f}});
  Array g = Array({{0.f}});
  Array b = Array({{0.f}});

  glm::vec3 target(1.f, 1.f, 1.f);

  // distance ≈ sqrt(3) ≈ 1.732
  Array mask = color_match_mask(r, g, b, target, 1.f);

  EXPECT_EQ(mask(0, 0), 0.f);
}
