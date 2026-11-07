#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(RecastCliff, CliffMaskFlatTerrainIsZero)
{
  Array z = Array(glm::ivec2(64, 64), 0.5f);
  Array cliff_mask;

  gpu::recast_cliff(z, 0.1f, 4, 0.2f, 2.f, 50, &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, 64);
  EXPECT_EQ(cliff_mask.shape.y, 64);
  EXPECT_NEAR(cliff_mask.max(), 0.f, 1e-6);
  EXPECT_NEAR(cliff_mask.min(), 0.f, 1e-6);
}

TEST(RecastCliff, CliffMaskSteepTerrainIsPositive)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float gain = 2.f;

  gpu::recast_cliff(z_copy, talus, ir, amplitude, gain, 50, &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);
  EXPECT_GT(cliff_mask.max(), 0.f);
}

TEST(RecastCliff, CliffMaskWithInputMask)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array mask = Array(shape, 0.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x / 2; ++i)
      mask(i, j) = 1.f;

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float gain = 2.f;

  gpu::recast_cliff(z_copy, talus, ir, amplitude, &mask, gain, 50, &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);

  // right half should be zero due to mask
  for (int j = 0; j < shape.y; ++j)
    for (int i = shape.x / 2; i < shape.x; ++i)
      EXPECT_NEAR(cliff_mask(i, j), 0.f, 1e-6);
}

TEST(RecastCliffDirectional, CliffMaskPositive)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float angle = 45.f;
  float gain = 2.f;

  gpu::recast_cliff_directional(z_copy,
                                talus,
                                ir,
                                amplitude,
                                angle,
                                gain,
                                50,
                                &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);
  EXPECT_GT(cliff_mask.max(), 0.f);
}

TEST(RecastCliffDirectional, CliffMaskWithInputMask)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array mask = Array(shape, 0.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x / 2; ++i)
      mask(i, j) = 1.f;

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float angle = 45.f;
  float gain = 2.f;

  gpu::recast_cliff_directional(z_copy,
                                talus,
                                ir,
                                amplitude,
                                angle,
                                &mask,
                                gain,
                                50,
                                &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);

  // right half should be zero due to mask
  for (int j = 0; j < shape.y; ++j)
    for (int i = shape.x / 2; i < shape.x; ++i)
      EXPECT_NEAR(cliff_mask(i, j), 0.f, 1e-6);
}

TEST(RecastCliffDirectional, CliffMaskWithVariableAngle)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z);

  Array angle_field = noise_fbm(NoiseType::PERLIN,
                                shape,
                                glm::vec2(1.f, 1.f),
                                10);
  remap(angle_field, 0.f, 360.f);

  Array z_copy = z;
  Array cliff_mask;
  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float gain = 2.f;

  gpu::recast_cliff_directional(z_copy,
                                talus,
                                ir,
                                amplitude,
                                angle_field,
                                gain,
                                50,
                                &cliff_mask);

  EXPECT_EQ(cliff_mask.shape.x, shape.x);
  EXPECT_EQ(cliff_mask.shape.y, shape.y);
  EXPECT_GE(cliff_mask.min(), 0.f);
  EXPECT_GT(cliff_mask.max(), 0.f);
}

TEST(RecastCliff, TiledDistributedExecution)
{
  const glm::ivec2 shape = {128, 128};
  const glm::ivec2 tile_shape = {32, 32};
  const int        halo = 8;

  Array z_arr = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z_arr);

  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float gain = 2.f;
  int   iterations = 40;
  int   sync_interval = 10;

  for (auto [mode, storage_mode] :
       {std::pair{ForEachMode::VA_SEQUENTIAL, StorageMode::VA_DISK_LRU_MIN},
        std::pair{ForEachMode::VA_DISTRIBUTED, StorageMode::VA_DISK_LRU},
        std::pair{ForEachMode::VA_DISTRIBUTED, StorageMode::VA_DISK_LRU_MIN}})
  {
    ComputeMode cm{.mode = mode, .trim_storage = true};

    VirtualArray z(shape, tile_shape, halo, storage_mode);
    z.from_array(z_arr, cm);

    for_each_tile(
        z,
        [&](Array &tile, const TileRegion &region, ICommunicator &comm)
        {
          if (mode == ForEachMode::VA_SEQUENTIAL)
          {
            EXPECT_LE(z.live_tile_count(), 2u);
          }

          gpu::recast_cliff(tile,
                            talus,
                            ir,
                            amplitude,
                            gain,
                            iterations,
                            nullptr,
                            sync_interval,
                            region,
                            comm);
        },
        cm);

    Array res = z.to_array(cm);

    EXPECT_GT(res.sum(), 0.f);
    EXPECT_FALSE(assert_almost_equal(res, z_arr, 1e-4f));

    if (mode == ForEachMode::VA_SEQUENTIAL)
    {
      EXPECT_LE(z.live_tile_count(), 2u);
    }
  }
}

TEST(RecastCliffDirectional, TiledDistributedExecution)
{
  const glm::ivec2 shape = {128, 128};
  const glm::ivec2 tile_shape = {32, 32};
  const int        halo = 8;

  Array z_arr = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(2.f, 2.f), 42);
  remap(z_arr);

  float talus = 1.f / shape.x;
  int   ir = 4;
  float amplitude = 0.1f;
  float angle = 45.f;
  float gain = 2.f;
  int   iterations = 40;
  int   sync_interval = 10;

  for (auto [mode, storage_mode] :
       {std::pair{ForEachMode::VA_SEQUENTIAL, StorageMode::VA_DISK_LRU_MIN},
        std::pair{ForEachMode::VA_DISTRIBUTED, StorageMode::VA_DISK_LRU},
        std::pair{ForEachMode::VA_DISTRIBUTED, StorageMode::VA_DISK_LRU_MIN}})
  {
    ComputeMode cm{.mode = mode, .trim_storage = true};

    VirtualArray z(shape, tile_shape, halo, storage_mode);
    z.from_array(z_arr, cm);

    for_each_tile(
        z,
        [&](Array &tile, const TileRegion &region, ICommunicator &comm)
        {
          if (mode == ForEachMode::VA_SEQUENTIAL)
          {
            EXPECT_LE(z.live_tile_count(), 2u);
          }

          gpu::recast_cliff_directional(tile,
                                        talus,
                                        ir,
                                        amplitude,
                                        angle,
                                        gain,
                                        iterations,
                                        nullptr,
                                        sync_interval,
                                        region,
                                        comm);
        },
        cm);

    Array res = z.to_array(cm);

    EXPECT_GT(res.sum(), 0.f);
    EXPECT_FALSE(assert_almost_equal(res, z_arr, 1e-4f));

    if (mode == ForEachMode::VA_SEQUENTIAL)
    {
      EXPECT_LE(z.live_tile_count(), 2u);
    }
  }
}
