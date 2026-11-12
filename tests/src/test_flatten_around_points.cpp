#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(FlattenAroundPoints, SinglePoint)
{
  hmap::Array z({65, 65}, 0.f);

  hmap::flatten_around_points(z, {0.5f}, {0.5f}, {1.f}, {0.2f});

  // exact reference at the point, untouched far away
  EXPECT_NEAR(z(32, 32), 1.f, 1e-5f);
  EXPECT_NEAR(z(0, 0), 0.f, 1e-6f);
}

TEST(FlattenAroundPoints, LargerRadiusDominates)
{
  hmap::Array z({65, 65}, 0.f);

  // cell at the middle is covered by both points; the large one dominates
  hmap::flatten_around_points(z,
                              {0.4f, 0.6f},
                              {0.5f, 0.5f},
                              {1.f, 2.f},
                              {0.1f, 0.4f});

  EXPECT_GT(z(32, 32), 1.5f);
}
