#include <cmath>
#include <vector>

#include "highmap/geometry/path.hpp"

#include <gtest/gtest.h>

using namespace hmap;

static constexpr float eps = 1e-5f;

TEST(PathScale, ScaleUniformDefaultCenter)
{
  Path path(std::vector<Point>{{0.f, 0.f, 1.f}, {1.f, 1.f, 2.f}});
  path.set_closed(false);

  Path scaled = scale(path, 0.8f);

  EXPECT_FALSE(scaled.is_closed());
  ASSERT_EQ(scaled.size(), 2u);
  EXPECT_NEAR(scaled[0].x, 0.1f, eps);
  EXPECT_NEAR(scaled[0].y, 0.1f, eps);
  EXPECT_NEAR(scaled[0].v, 1.f, eps);

  EXPECT_NEAR(scaled[1].x, 0.9f, eps);
  EXPECT_NEAR(scaled[1].y, 0.9f, eps);
  EXPECT_NEAR(scaled[1].v, 2.f, eps);
}

TEST(PathScale, ScaleClosedPathPreservesClosure)
{
  Path path(std::vector<Point>{{0.f, 0.f, 0.f},
                               {1.f, 0.f, 0.f},
                               {1.f, 1.f, 0.f},
                               {0.f, 1.f, 0.f}});
  path.set_closed(true);

  Path scaled = scale(path, 0.5f);

  EXPECT_TRUE(scaled.is_closed());
  ASSERT_EQ(scaled.size(), 4u);
  EXPECT_NEAR(scaled[0].x, 0.25f, eps);
  EXPECT_NEAR(scaled[0].y, 0.25f, eps);
  EXPECT_NEAR(scaled[2].x, 0.75f, eps);
  EXPECT_NEAR(scaled[2].y, 0.75f, eps);
}

TEST(PathScale, ScaleNonUniformCustomCenter)
{
  Path path(std::vector<Point>{{1.f, 2.f, 3.f}, {3.f, 4.f, 5.f}});

  Path scaled = scale(path, glm::vec2(2.f, 0.5f), glm::vec2(1.f, 2.f));

  ASSERT_EQ(scaled.size(), 2u);
  // Point 0: at center (1, 2) -> (1, 2)
  EXPECT_NEAR(scaled[0].x, 1.f, eps);
  EXPECT_NEAR(scaled[0].y, 2.f, eps);
  EXPECT_NEAR(scaled[0].v, 3.f, eps);

  // Point 1: (3, 4) -> x = 1 + 2 * (3 - 1) = 5, y = 2 + 0.5 * (4 - 2) = 3
  EXPECT_NEAR(scaled[1].x, 5.f, eps);
  EXPECT_NEAR(scaled[1].y, 3.f, eps);
  EXPECT_NEAR(scaled[1].v, 5.f, eps);
}

TEST(PathScale, ScaleExpandBeyondDomain)
{
  Path path(std::vector<Point>{{0.2f, 0.2f, 0.f}, {0.8f, 0.8f, 0.f}});

  // Scale 2.0 with default center (0.5, 0.5)
  // p0: 0.5 + 2.0 * (0.2 - 0.5) = -0.1
  // p1: 0.5 + 2.0 * (0.8 - 0.5) = 1.1
  Path scaled = scale(path, 2.0f);

  ASSERT_EQ(scaled.size(), 2u);
  EXPECT_NEAR(scaled[0].x, -0.1f, eps);
  EXPECT_NEAR(scaled[0].y, -0.1f, eps);
  EXPECT_NEAR(scaled[1].x, 1.1f, eps);
  EXPECT_NEAR(scaled[1].y, 1.1f, eps);
}

TEST(PathScale, ScalePreservesPointValues)
{
  Path path(std::vector<Point>{
      {0.2f, 0.3f, 10.0f},
      {0.4f, 0.6f, 20.0f},
      {0.7f, 0.8f, 30.0f},
  });

  Path scaled = scale(path, 0.5f);

  ASSERT_EQ(scaled.size(), 3u);
  EXPECT_FLOAT_EQ(scaled[0].v, 10.0f);
  EXPECT_FLOAT_EQ(scaled[1].v, 20.0f);
  EXPECT_FLOAT_EQ(scaled[2].v, 30.0f);
}
