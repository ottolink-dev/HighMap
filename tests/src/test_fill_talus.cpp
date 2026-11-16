#include "highmap/filters.hpp"
#include "highmap/internal/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(FillTalus, DeterministicWithSameSeed)
{
  Array z1 = Array({{1, 1, 1}, {1, 5, 1}, {1, 1, 1}});

  Array z2 = z1;

  Array talus = Array(z1.shape, 1.f);

  fill_talus(z1, talus, 123, 1, 0.2f, nullptr);
  fill_talus(z2, talus, 123, 1, 0.2f, nullptr);

  EXPECT_TRUE(assert_almost_equal(z1, z2));
}

TEST(FillTalus, DifferentSeedsProduceDifferentResults)
{
  Array z1 = Array({{1, 1, 1}, {1, 5, 1}, {1, 1, 1}});

  Array z2 = z1;

  Array talus = Array(z1.shape, 1.f);

  fill_talus(z1, talus, 1, 1, 0.5f, nullptr);
  fill_talus(z2, talus, 2, 1, 0.5f, nullptr);

  // not strict equality → just ensure difference exists
  bool different = false;

  for (int i = 0; i < z1.shape.x; ++i)
    for (int j = 0; j < z1.shape.y; ++j)
      if (std::abs(z1(i, j) - z2(i, j)) > 1e-5f) different = true;

  EXPECT_TRUE(different);
}

TEST(FillTalus, SeedMaskRestrictsPropagation)
{
  Array z = Array({{1, 1, 1}, {1, 10, 1}, {1, 1, 1}});

  Array talus = Array(z.shape, 1.f);

  Array mask = Array(z.shape, 0.f);
  mask(1, 1) = 1.f; // only center is seed

  fill_talus(z, talus, 42, 1, 0.f, &mask);

  // neighbors should increase, corners may remain lower
  EXPECT_GT(z(0, 1), 1.f);
  EXPECT_GT(z(1, 0), 1.f);
}

TEST(FillTalus, RespectsLocalSlopeConstraint)
{
  Array z = Array({{0, 0, 0}, {0, 10, 0}, {0, 0, 0}});

  Array talus = Array(z.shape, 1.f);

  fill_talus(z, talus, 42, 1, 0.f, nullptr);

  for (int i = 0; i < z.shape.x; ++i)
    for (int j = 0; j < z.shape.y; ++j)
      for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj)
        {
          int ni = i + di;
          int nj = j + dj;

          if (ni < 0 || nj < 0 || ni >= z.shape.x || nj >= z.shape.y) continue;

          float dist = std::hypot(di, dj);
          if (dist == 0.f) continue;

          float slope = z(i, j) - z(ni, nj);

          EXPECT_LE(slope, dist * talus(i, j) + 1e-4f);
        }
}
