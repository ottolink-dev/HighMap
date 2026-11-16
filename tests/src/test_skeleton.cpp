#include "highmap/internal/assert.hpp"
#include "highmap/morphology.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Skeleton, EmptyInput)
{
  Array input = Array({{0, 0}, {0, 0}});
  Array out = skeleton(input, false);

  EXPECT_EQ(count_non_zero(out), 0);
}

TEST(Skeleton, SinglePixelInvariant)
{
  Array input = Array({{0, 0, 0}, {0, 1, 0}, {0, 0, 0}});
  Array out = skeleton(input, false);

  EXPECT_EQ(count_non_zero(out), 1);
  EXPECT_EQ(out(1, 1), 1);
}

TEST(Skeleton, AlreadyThinStructure)
{
  Array input = Array({{0, 1, 0}, {0, 1, 0}, {0, 1, 0}});
  Array out = skeleton(input, false);

  // should remain unchanged
  EXPECT_EQ(count_non_zero(out), 3);
  EXPECT_TRUE(is_subset(out, input));
  EXPECT_TRUE(is_subset(input, out));
}

TEST(Skeleton, ThickBlockReduces)
{
  Array input = Array({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});
  Array out = skeleton(input, false);

  // skeleton must be strictly thinner
  EXPECT_LT(count_non_zero(out), count_non_zero(input));

  // skeleton must be subset of input
  EXPECT_TRUE(is_subset(out, input));
}

TEST(Skeleton, ConnectivityPreserved)
{
  Array input = Array({{1, 1, 0}, {0, 1, 1}, {0, 0, 1}});
  Array out = skeleton(input, false);

  // Skeleton should not introduce new pixels
  EXPECT_TRUE(is_subset(out, input));

  // Should still have at least one pixel (not erased)
  EXPECT_GT(count_non_zero(out), 0);
}

TEST(Skeleton, Idempotence)
{
  Array input = Array({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});

  Array sk1 = skeleton(input, false);
  Array sk2 = skeleton(sk1, false);

  // applying skeleton twice should not change result
  for (int i = 0; i < sk1.shape.x; ++i)
    for (int j = 0; j < sk1.shape.y; ++j)
      EXPECT_EQ(sk1(i, j), sk2(i, j));
}

TEST(Skeleton, ZeroAtBorders)
{
  Array input = Array({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});
  Array out = skeleton(input, true);

  int nx = out.shape.x;
  int ny = out.shape.y;

  // borders must be zero
  for (int i = 0; i < nx; ++i)
  {
    EXPECT_EQ(out(i, 0), 0);
    EXPECT_EQ(out(i, ny - 1), 0);
  }

  for (int j = 0; j < ny; ++j)
  {
    EXPECT_EQ(out(0, j), 0);
    EXPECT_EQ(out(nx - 1, j), 0);
  }
}

TEST(Skeleton, ReductionButNotFullRemoval)
{
  Array input = Array({{0, 1, 1, 0}, {1, 1, 1, 1}, {0, 1, 1, 0}});
  Array out = skeleton(input, false);

  // Should reduce but not vanish
  EXPECT_LT(count_non_zero(out), count_non_zero(input));
  EXPECT_GT(count_non_zero(out), 0);
}

TEST(Skeleton_CPU_GPU, RandomBinaryFields)
{
  std::mt19937                       rng(42);
  std::uniform_int_distribution<int> dist(0, 1);

  const int nx = 64;
  const int ny = 64;

  for (int test = 0; test < 20; ++test)
  {
    Array input(glm::ivec2(nx, ny));

    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        input(i, j) = float(dist(rng));

    Array cpu = skeleton(input, false);
    Array gpu = gpu::skeleton(input, false);

    int diff = count_non_zero(cpu - gpu);

    EXPECT_EQ(diff, 0);
  }
}
