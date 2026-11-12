#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(FlattenAroundPointsTerrainAware, FlatTerrainPlateau)
{
  hmap::Array z({65, 65}, 0.f);

  hmap::flatten_around_points_terrain_aware(z,
                                            {0.5f},
                                            {0.5f},
                                            {1.f},
                                            0.2f,
                                            1.f,
                                            2.f,
                                            -1.f,
                                            0.f);

  EXPECT_NEAR(z(32, 32), 1.f, 1e-5f);
  EXPECT_NEAR(z(0, 0), 0.f, 1e-6f);
}

TEST(FlattenAroundPointsTerrainAware, FollowsValley)
{
  // valley along y (low at x = 0.5), walls rising with |x - 0.5|
  const int   n = 65;
  hmap::Array z({n, n}, 0.f);

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      z(i, j) = 4.f * std::abs(float(i) / float(n - 1) - 0.5f);

  hmap::Array mask;
  hmap::flatten_around_points_terrain_aware(z,
                                            {0.5f},
                                            {0.5f},
                                            {0.f},
                                            0.3f,
                                            5.f,
                                            20.f,
                                            -1.f,
                                            0.f,
                                            1.f,
                                            hmap::RadialProfile::RP_SMOOTHSTEP,
                                            1.f,
                                            &mask);

  // reaches further along the valley than across it
  EXPECT_GT(mask(32, 32 + 12), mask(32 + 12, 32));
}

TEST(FlattenAroundPointsTerrainAware, OutsideBboxSeedIgnored)
{
  hmap::Array z({33, 33}, 2.f);

  hmap::flatten_around_points_terrain_aware(z, {2.f}, {2.f}, {0.f}, 0.2f);

  EXPECT_NEAR(z(16, 16), 2.f, 1e-6f);
}
