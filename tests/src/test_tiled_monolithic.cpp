#include "highmap/dbg/assert.hpp"
#include "highmap/erosion.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// --- Depression filling tests

TEST(DepressionFilling_Communicator, TiledMatchesMonolithic)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  // monolithic reference
  Array ref = input;
  depression_filling(ref, 200, 1e-4f);

  for (auto mode : {ForEachMode::VA_SEQUENTIAL, ForEachMode::VA_DISTRIBUTED})
  {
    VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
    ComputeMode  cm{.mode = mode};

    va.from_array(input, cm);

    for_each_tile(
        va,
        [](Array &tile, const TileRegion &region, ICommunicator &comm)
        {
          depression_filling(tile,
                             200,
                             1e-4f,
                             true,
                             true,
                             true,
                             true,
                             region,
                             comm);
        },
        cm);

    Array result = va.to_array(cm);
    Array diff = abs(result - ref);

    // both fill depressions and produce elevations >= input
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        EXPECT_GE(result(i, j), input(i, j) - 1e-4f);
      }

    EXPECT_LT(diff.mean(), 0.02f)
        << "DepressionFilling mean diff exceeded for mode=" << to_string(mode);
    EXPECT_LT(diff.max(), 0.1f)
        << "DepressionFilling max diff exceeded for mode=" << to_string(mode);
  }
}

// --- Priority-Flood depression filling tests

TEST(DepressionFillingPriorityFlood_Communicator, TiledMatchesMonolithic)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  // monolithic reference
  Array ref = input;
  depression_filling_priority_flood(ref, false);

  for (auto mode : {ForEachMode::VA_SEQUENTIAL, ForEachMode::VA_DISTRIBUTED})
  {
    VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
    ComputeMode  cm{.mode = mode};

    va.from_array(input, cm);

    for_each_tile(
        va,
        [](Array &tile, const TileRegion &region, ICommunicator &comm)
        {
          depression_filling_priority_flood(tile,
                                            false,
                                            true,
                                            true,
                                            true,
                                            true,
                                            region,
                                            comm);
        },
        cm);

    Array result = va.to_array(cm);
    Array diff = abs(result - ref);

    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        EXPECT_GE(result(i, j), input(i, j) - 1e-4f);
      }

    EXPECT_LT(diff.mean(), 0.01f)
        << "DepressionFillingPriorityFlood mean diff exceeded for mode="
        << to_string(mode);
    EXPECT_LT(diff.max(), 0.05f)
        << "DepressionFillingPriorityFlood max diff exceeded for mode="
        << to_string(mode);
  }
}

// --- Flooding lake system tests

TEST(FloodingLakeSystem_VirtualArray, SingleTileMatchesArray)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{64, 64};
  const int        halo = 0;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};

  va.from_array(input, cm);

  VirtualArray va_water;
  va_water.copy_from(va, cm);

  hmap::for_each_tile(
      va_water,
      [](Array &tile, const TileRegion &region, ICommunicator &comm)
      { tile = flooding_lake_system(tile, 0.f, region, comm); },
      cm);

  Array expected = flooding_lake_system(input, 0.f);
  Array result = va_water.to_array(cm);

  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-5f));
}

TEST(FloodingLakeSystem_VirtualArray, MultiTileWithSurfaceThreshold)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};

  va.from_array(input, cm);

  VirtualArray va_water;
  va_water.copy_from(va, cm);

  hmap::for_each_tile(
      va_water,
      [](Array &tile, const TileRegion &region, ICommunicator &comm)
      { tile = flooding_lake_system(tile, 10.f, region, comm); },
      cm);

  Array result = va_water.to_array(cm);

  EXPECT_GE(result.min(), 0.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      EXPECT_GE(result(i, j), 0.f);
}

// --- Flooding from point tests

TEST(FloodingFromPoint_Communicator, TiledMatchesMonolithicSingleSeed)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {3.f, 3.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 123);

  remap(input, 0.f, 10.f);

  int   seed_x = 20;
  int   seed_y = 25;
  float depth_min = 2.0f;

  Array expected = flooding_from_point(input, seed_x, seed_y, depth_min);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};

  va.from_array(input, cm);

  VirtualArray va_water;
  va_water.copy_from(va, cm);

  hmap::for_each_tile(
      va_water,
      [seed_x, seed_y, depth_min](Array            &tile,
                                  const TileRegion &region,
                                  ICommunicator    &comm)
      {
        tile = flooding_from_point(tile,
                                   seed_x,
                                   seed_y,
                                   depth_min,
                                   region,
                                   comm);
      },
      cm);

  Array result = va_water.to_array(cm);

  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-4f));
}

TEST(FloodingFromPoint_Communicator, TiledMatchesMonolithicBoundarySeed)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {3.f, 3.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 456);

  remap(input, 0.f, 10.f);

  // seed placed directly on the tile boundary (x = 31, y = 32)
  int   seed_x = 31;
  int   seed_y = 32;
  float depth_min = 3.0f;

  Array expected = flooding_from_point(input, seed_x, seed_y, depth_min);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};

  va.from_array(input, cm);

  VirtualArray va_water;
  va_water.copy_from(va, cm);

  hmap::for_each_tile(
      va_water,
      [seed_x, seed_y, depth_min](Array            &tile,
                                  const TileRegion &region,
                                  ICommunicator    &comm)
      {
        tile = flooding_from_point(tile,
                                   seed_x,
                                   seed_y,
                                   depth_min,
                                   region,
                                   comm);
      },
      cm);

  Array result = va_water.to_array(cm);

  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-4f));
}

TEST(FloodingFromPoint_Communicator, TiledMatchesMonolithicMultiSeed)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {3.f, 3.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 789);

  remap(input, 0.f, 10.f);

  std::vector<int> seeds_x = {10, 31, 48};
  std::vector<int> seeds_y = {15, 32, 50};
  float            depth_min = 2.5f;

  Array expected = flooding_from_point(input, seeds_x, seeds_y, depth_min);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};

  va.from_array(input, cm);

  VirtualArray va_water;
  va_water.copy_from(va, cm);

  hmap::for_each_tile(
      va_water,
      [&seeds_x, &seeds_y, depth_min](Array            &tile,
                                      const TileRegion &region,
                                      ICommunicator    &comm)
      {
        tile = flooding_from_point(tile,
                                   seeds_x,
                                   seeds_y,
                                   depth_min,
                                   region,
                                   comm);
      },
      cm);

  Array result = va_water.to_array(cm);

  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-4f));
}

// --- Flooding from boundaries tests

TEST(FloodingFromBoundaries_Communicator, TiledMatchesMonolithicAllBoundaries)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {3.f, 3.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 101);

  remap(input, 0.f, 10.f);

  float zref = 6.0f;
  Array expected = flooding_from_boundaries(input,
                                            zref,
                                            true,
                                            true,
                                            true,
                                            true);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};

  va.from_array(input, cm);

  VirtualArray va_water;
  va_water.copy_from(va, cm);

  hmap::for_each_tile(
      va_water,
      [zref](Array &tile, const TileRegion &region, ICommunicator &comm)
      {
        tile = flooding_from_boundaries(tile,
                                        zref,
                                        true,
                                        true,
                                        true,
                                        true,
                                        region,
                                        comm);
      },
      cm);

  Array result = va_water.to_array(cm);

  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-4f));
}

TEST(FloodingFromBoundaries_Communicator,
     TiledMatchesMonolithicSelectiveBoundaries)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {3.f, 3.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 202);

  remap(input, 0.f, 10.f);

  // flood only from west and south boundaries
  float zref = 5.5f;
  Array expected = flooding_from_boundaries(input,
                                            zref,
                                            false,
                                            true,
                                            false,
                                            true);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};

  va.from_array(input, cm);

  VirtualArray va_water;
  va_water.copy_from(va, cm);

  hmap::for_each_tile(
      va_water,
      [zref](Array &tile, const TileRegion &region, ICommunicator &comm)
      {
        tile = flooding_from_boundaries(tile,
                                        zref,
                                        false,
                                        true,
                                        false,
                                        true,
                                        region,
                                        comm);
      },
      cm);

  Array result = va_water.to_array(cm);

  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-4f));
}
