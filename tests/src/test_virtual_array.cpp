#include "highmap/dbg/assert.hpp"
#include "highmap/primitives.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(VirtualArrayTest, ConstructionAndProperties)
{
  glm::ivec2 shape{100, 80};
  glm::vec4  bbox{0.f, 10.f, 0.f, 8.f};
  glm::ivec2 tile_shape{32, 32};
  int        halo = 4;

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);

  EXPECT_EQ(va.shape, shape);
  EXPECT_EQ(va.bbox, bbox);
  EXPECT_EQ(va.tile_shape, tile_shape);
  EXPECT_EQ(va.halo, halo);

  glm::ivec2 grid = va.get_max_tiles();
  EXPECT_EQ(grid.x, 4); // ceil(100/32) = 4
  EXPECT_EQ(grid.y, 3); // ceil(80/32) = 3
  EXPECT_EQ(va.get_ntiles(), 12);
}

TEST(VirtualArrayTest, SingleTileGeometry)
{
  glm::ivec2   shape{16, 16};
  glm::ivec2   tile_shape{32, 32};
  int          halo = 2;
  VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);

  EXPECT_EQ(va.get_ntiles(), 1);
  TileRegion region = va.tile_region_from_tile_coords(0, 0);
  EXPECT_EQ(region.shape.x, 16);
  EXPECT_EQ(region.shape.y, 16);
  // halos at outer boundaries are 0
  EXPECT_EQ(region.halo.x, 0);
  EXPECT_EQ(region.halo.y, 0);
  EXPECT_EQ(region.halo.z, 0);
  EXPECT_EQ(region.halo.w, 0);
}

TEST(VirtualArrayTest, TileCoordinateMapping)
{
  glm::ivec2   shape{64, 64};
  glm::ivec2   tile_shape{32, 32};
  int          halo = 4;
  VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);

  TileRegion r00 = va.tile_region_from_tile_coords(0, 0);
  EXPECT_EQ(r00.halo.x, 0);
  EXPECT_EQ(r00.halo.y, 4);
  EXPECT_EQ(r00.halo.z, 0);
  EXPECT_EQ(r00.halo.w, 4);
  EXPECT_EQ(r00.shape.x, 36);
  EXPECT_EQ(r00.shape.y, 36);

  TileRegion r11 = va.tile_region_from_tile_coords(1, 1);
  EXPECT_EQ(r11.halo.x, 4);
  EXPECT_EQ(r11.halo.y, 0);
  EXPECT_EQ(r11.halo.z, 4);
  EXPECT_EQ(r11.halo.w, 0);

  TileRegion r_at = va.tile_region_from_global_index(40, 10);
  EXPECT_EQ(r_at.key.tx, 1);
  EXPECT_EQ(r_at.key.ty, 0);
}

TEST(VirtualArrayTest, CellGetSet)
{
  glm::ivec2   shape{32, 32};
  glm::ivec2   tile_shape{16, 16};
  int          halo = 2;
  VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);

  va.set(5, 5, 42.0f);
  va.set(20, 25, 99.0f);

  EXPECT_FLOAT_EQ(va.get(5, 5), 42.0f);
  EXPECT_FLOAT_EQ(va.get(20, 25), 99.0f);
  EXPECT_FLOAT_EQ(va.get(0, 0), 0.0f);
}

TEST(VirtualArrayTest, FillAndReductions)
{
  glm::ivec2   shape{40, 40};
  glm::ivec2   tile_shape{16, 16};
  VirtualArray va(shape, tile_shape, 2, StorageMode::VA_RAM);

  ComputeMode cm{.mode = ForEachMode::VA_SEQUENTIAL};
  va.fill(3.5f, cm);

  EXPECT_FLOAT_EQ(va.min(cm), 3.5f);
  EXPECT_FLOAT_EQ(va.max(cm), 3.5f);
  EXPECT_FLOAT_EQ(va.mean(cm), 3.5f);
  EXPECT_FLOAT_EQ(va.sum(cm), 3.5f * 40 * 40);

  glm::vec2 rng = va.range(cm);
  EXPECT_FLOAT_EQ(rng.x, 3.5f);
  EXPECT_FLOAT_EQ(rng.y, 3.5f);
}

TEST(VirtualArrayTest, ArrayRoundTripSequentialAndDistributed)
{
  glm::ivec2 shape{48, 36};
  Array      src(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      src(i, j) = float(i + 2 * j);

  for (auto mode : {ForEachMode::VA_SEQUENTIAL, ForEachMode::VA_DISTRIBUTED})
  {
    ComputeMode  cm{.mode = mode};
    VirtualArray va(shape, {16, 16}, 4, StorageMode::VA_RAM);

    va.from_array(src, cm);
    Array dst = va.to_array(cm);

    EXPECT_EQ(dst.shape, shape);
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        EXPECT_FLOAT_EQ(dst(i, j), src(i, j))
            << "Mismatch at (" << i << ", " << j << ") with mode "
            << to_string(mode);
      }
  }
}

TEST(VirtualArrayTest, CloneAndCopyFrom)
{
  glm::ivec2   shape{32, 32};
  VirtualArray va(shape, {16, 16}, 2, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};
  va.fill(12.34f, cm);

  auto cloned = va.clone(cm, /*deep_copy=*/true);
  EXPECT_EQ(cloned->shape, va.shape);
  EXPECT_NEAR(cloned->mean(cm), 12.34f, 1e-4f);

  VirtualArray va_copy(shape, {16, 16}, 2, StorageMode::VA_RAM);
  va_copy.copy_from(va, cm, /*copy_src_data=*/true);
  EXPECT_NEAR(va_copy.mean(cm), 12.34f, 1e-4f);
}

TEST(VirtualArrayTest, RemapAndInverse)
{
  glm::ivec2   shape{32, 32};
  VirtualArray va(shape, {16, 16}, 2, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};

  // set gradient
  for_each_tile(
      va,
      [](Array &tile, const TileRegion &region)
      {
        for (int j = 0; j < region.shape.y; ++j)
          for (int i = 0; i < region.shape.x; ++i)
            tile(i, j) = float(region.key.tx + region.key.ty);
      },
      cm);

  va.remap(0.f, 1.f, cm);
  EXPECT_NEAR(va.min(cm), 0.f, 1e-5f);
  EXPECT_NEAR(va.max(cm), 1.f, 1e-5f);

  va.inverse(cm);
  EXPECT_NEAR(va.min(cm), 0.f, 1e-5f);
  EXPECT_NEAR(va.max(cm), 1.f, 1e-5f);
}

TEST(VirtualArrayTest, StorageModesPersistence)
{
  glm::ivec2  shape{32, 32};
  glm::ivec2  tile_shape{16, 16};
  ComputeMode cm{.mode = ForEachMode::VA_SEQUENTIAL};

  for (auto mode : {StorageMode::VA_RAM,
                    StorageMode::VA_DISK_LRU,
                    StorageMode::VA_DISK_LRU_MIN,
                    StorageMode::VA_DISK_SEQUENTIAL})
  {
    VirtualArray va(shape, tile_shape, 2, mode);
    va.fill(7.0f, cm);

    EXPECT_FLOAT_EQ(va.min(cm), 7.0f)
        << "Failed for storage " << to_string(mode);
    EXPECT_FLOAT_EQ(va.max(cm), 7.0f)
        << "Failed for storage " << to_string(mode);
    EXPECT_FLOAT_EQ(va.mean(cm), 7.0f)
        << "Failed for storage " << to_string(mode);

    va.trim_storage();
    // After trim, data should still be loadable (except sequential which loads
    // per tile)
    if (mode != StorageMode::VA_RAM)
    {
      EXPECT_FLOAT_EQ(va.mean(cm), 7.0f)
          << "Failed after trim for storage " << to_string(mode);
    }
  }
}

TEST(VirtualArrayTest, MemoryFootprintTracking)
{
  glm::ivec2  shape{64, 64};
  glm::ivec2  tile_shape{16, 16}; // 4x4 = 16 tiles total
  ComputeMode cm{.mode = ForEachMode::VA_SEQUENTIAL};

  // 1. RAM storage holds all 16 tiles
  {
    VirtualArray va(shape, tile_shape, 2, StorageMode::VA_RAM);
    EXPECT_EQ(va.live_tile_count(), 0);
    EXPECT_EQ(va.live_memory_bytes(), 0);

    va.fill(1.0f, cm);
    EXPECT_EQ(va.live_tile_count(), 16);
    EXPECT_GT(va.live_memory_bytes(), 0);
  }

  // 2. DISK_LRU_MIN holds at most 2 tiles
  {
    VirtualArray va(shape, tile_shape, 2, StorageMode::VA_DISK_LRU_MIN);
    va.fill(2.0f, cm);
    EXPECT_LE(va.live_tile_count(), 2);

    va.trim_storage();
    EXPECT_EQ(va.live_tile_count(), 0);
    EXPECT_EQ(va.live_memory_bytes(), 0);
  }

  // 3. DISK_SEQUENTIAL holds at most 1 tile
  {
    VirtualArray va(shape, tile_shape, 2, StorageMode::VA_DISK_SEQUENTIAL);
    va.fill(3.0f, cm);
    EXPECT_LE(va.live_tile_count(), 1);

    va.trim_storage();
    EXPECT_EQ(va.live_tile_count(), 0);
    EXPECT_EQ(va.live_memory_bytes(), 0);
  }
}
