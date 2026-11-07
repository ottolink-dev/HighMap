#include "highmap/dbg/assert.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/math/array.hpp"
#include "highmap/morphology.hpp"
#include "highmap/primitives.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(DistanceTransform, ProducesCorrectDistances)
{
  glm::ivec2    shape = {128, 64};
  std::uint32_t seed = 42;

  // --- Generate random mask

  Array mask = white_sparse(shape, 0.f, 1.f, 0.05f, seed);

  // --- Brute force distance transform

  Array ref(shape);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float best = std::numeric_limits<float>::max();

      for (int jj = 0; jj < shape.y; ++jj)
        for (int ii = 0; ii < shape.x; ++ii)
        {
          if (mask(ii, jj) > 0.f)
          {
            float dx = float(i - ii);
            float dy = float(j - jj);
            float d2 = dx * dx + dy * dy;

            if (d2 < best) best = d2;
          }
        }

      ref(i, j) = std::sqrt(best);
    }

  // --- Fast EDT
  Array edt = distance_transform(mask);
  Array edt_jfa = gpu::distance_transform_jfa(mask);

  EXPECT_TRUE(assert_almost_equal(ref, edt, 1e-6f));
  EXPECT_TRUE(assert_almost_equal(ref, edt_jfa, 1e-3f));
}

TEST(DistanceTransformManhattan, ProducesCorrectDistances)
{
  glm::ivec2    shape = {128, 64};
  std::uint32_t seed = 42;

  // --- Generate random mask

  Array mask = white_sparse(shape, 0.f, 1.f, 0.05f, seed);

  // --- Brute force distance transform

  Array ref(shape);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float best = std::numeric_limits<float>::max();

      for (int jj = 0; jj < shape.y; ++jj)
        for (int ii = 0; ii < shape.x; ++ii)
        {
          if (mask(ii, jj) > 0.f)
          {
            float d = std::abs(i - ii) + std::abs(j - jj);
            best = std::min(best, d);
          }
        }

      ref(i, j) = best;
    }

  // --- Fast EDT

  Array edt = distance_transform_manhattan(mask);

  EXPECT_TRUE(assert_almost_equal(ref, edt, 1e-6f));
}

TEST(DistanceTransform, SquaredDistanceConsistency)
{
  glm::ivec2    shape = {128, 64};
  std::uint32_t seed = 42;

  Array mask = white_sparse(shape, 0.f, 1.f, 0.05f, seed);

  {
    Array edt = distance_transform(mask);
    Array edt2 = distance_transform(mask, true);
    EXPECT_TRUE(assert_almost_equal(edt * edt, edt2, 1e-6f));
  }

  {
    Array edt = distance_transform_manhattan(mask);
    Array edt2 = distance_transform_manhattan(mask, true);
    EXPECT_TRUE(assert_almost_equal(edt * edt, edt2, 1e-6f));
  }

  {
    Array edt = distance_transform_approx(mask);
    Array edt2 = distance_transform_approx(mask, true);
    EXPECT_TRUE(assert_almost_equal(edt * edt, edt2, 1e-6f));
  }
}

TEST(DistanceTransform, EmptyMask)
{
  glm::ivec2 shape = {128, 64};
  Array      mask(shape, 0.f);
  Array      edt = distance_transform(mask);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      ASSERT_TRUE(edt(i, j) >= shape.x + shape.y - 1.f)
          << "Expected large/infinite distance at (" << i << "," << j << ")";
    }
}

TEST(DistanceTransform, FullMask)
{
  glm::ivec2 shape = {128, 64};
  Array      mask(shape, 1.f);
  Array      edt = distance_transform(mask);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      ASSERT_NEAR(edt(i, j), 0.0f, 1e-6f)
          << "Expected zero distance at (" << i << "," << j << ")";
    }
}

TEST(DistanceTransform, SinglePoint)
{
  glm::ivec2 shape = {64, 64};

  Array mask(shape, 0.f);

  int cx = shape.x / 2;
  int cy = shape.y / 2;

  mask(cx, cy) = 1.f;

  Array edt = distance_transform(mask);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = float(i - cx);
      float dy = float(j - cy);
      float expected = std::sqrt(dx * dx + dy * dy);

      ASSERT_NEAR(edt(i, j), expected, 1e-4f)
          << "Mismatch at (" << i << "," << j << ")";
    }
}

TEST(DistanceTransform, HorizontalLine)
{
  glm::ivec2 shape = {64, 64};

  Array mask(shape);

  int line_j = shape.y / 2;

  for (int i = 0; i < shape.x; ++i)
    mask(i, line_j) = 1.f;

  Array edt = distance_transform(mask);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float expected = std::abs(float(j - line_j));

      ASSERT_NEAR(edt(i, j), expected, 1e-4f)
          << "Mismatch at (" << i << "," << j << ")";
    }
}

TEST(DistanceTransform, VerticalLine)
{
  glm::ivec2 shape = {64, 64};

  Array mask(shape);

  int line_i = shape.x / 2;

  for (int j = 0; j < shape.y; ++j)
    mask(line_i, j) = 1.f;

  Array edt = distance_transform(mask);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float expected = std::abs(float(i - line_i));

      ASSERT_NEAR(edt(i, j), expected, 1e-4f)
          << "Mismatch at (" << i << "," << j << ")";
    }
}

TEST(DistanceTransform, ConsistencyWithDistanceTransformWithClosest)
{
  std::uint32_t   seed = 0;
  glm::ivec2      shape = {64, 32};
  Mat<glm::ivec2> closest(shape); // not used

  for (int test = 0; test < 20; ++test)
  {
    Array mask = white_sparse(shape, 0.f, 1.f, 0.05f, seed++);
    Array a = distance_transform(mask);
    Array b = distance_transform_with_closest(mask, closest);

    int diff = count_non_zero(a - b);

    EXPECT_EQ(diff, 0);
  }
}

TEST(DistanceTransform, ApproxMatchesExactWithinTolerance)
{
  glm::ivec2    shape = {128, 64};
  std::uint32_t seed = 42;

  Array mask = white_sparse(shape, 0.f, 1.f, 0.05f, seed);

  Array exact = distance_transform(mask);
  Array approx = distance_transform_approx(mask, false);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float e = exact(i, j);
      float a = approx(i, j);

      // allow approximation error
      ASSERT_NEAR(a, e, 3.0f) << "Mismatch at (" << i << "," << j << ")";
    }
}

TEST(DistanceTransform, TiledWithCommunicator)
{
  const glm::ivec2 shape{128, 128};
  const glm::ivec2 tile_shape{64, 64};
  const int        halo = 8;

  // Test 1: Random sparse mask across 2x2 tiles
  {
    Array mask = white_sparse(shape, 0.f, 1.f, 0.05f, 42);
    Array ref = distance_transform(mask);

    VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);
    ComputeMode  cm_dist{.mode = ForEachMode::VA_DISTRIBUTED};
    va.from_array(mask, cm_dist);

    for_each_tile(
        va,
        [](Array &tile, const TileRegion &region, ICommunicator &comm)
        { tile = distance_transform(tile, false, region, comm); },
        cm_dist);

    Array res = va.to_array(cm_dist);
    Array diff = abs(res - ref);

    EXPECT_LT(diff.max(), 1e-4f);
  }

  // Test 2: Single seed in tile (0, 0) propagating into empty neighboring tiles
  {
    Array mask(shape, 0.f);
    mask(16, 16) = 1.f;

    Array ref = distance_transform(mask);

    VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);
    ComputeMode  cm_dist{.mode = ForEachMode::VA_DISTRIBUTED};
    va.from_array(mask, cm_dist);

    for_each_tile(
        va,
        [](Array &tile, const TileRegion &region, ICommunicator &comm)
        { tile = distance_transform(tile, false, region, comm); },
        cm_dist);

    Array res = va.to_array(cm_dist);
    Array diff = abs(res - ref);

    EXPECT_LT(diff.max(), 1e-4f);
  }

  // Test 3: Single seed in bottom-right tile (1, 1) propagating West and North
  // in VA_SEQUENTIAL mode with VA_DISK_LRU_MIN (2 live tiles)
  {
    Array mask(shape, 0.f);
    mask(100, 100) = 1.f;

    Array ref = distance_transform(mask);

    VirtualArray va(shape, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);
    ComputeMode  cm_seq{.mode = ForEachMode::VA_SEQUENTIAL,
                        .trim_storage = true};
    va.from_array(mask, cm_seq);

    for_each_tile(
        va,
        [](Array &tile, const TileRegion &region, ICommunicator &comm)
        { tile = distance_transform(tile, false, region, comm); },
        cm_seq);

    EXPECT_LE(va.live_tile_count(), 2u);

    Array res = va.to_array(cm_seq);
    Array diff = abs(res - ref);

    EXPECT_LT(diff.max(), 1e-4f);
  }
}
