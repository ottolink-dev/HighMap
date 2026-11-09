#include "highmap/dbg/assert.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
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

TEST(VirtualArrayTest, TileRegionInnerBbox)
{
  glm::ivec2   shape{64, 64};
  glm::vec4    bbox{0.f, 64.f, 0.f, 64.f};
  glm::ivec2   tile_shape{32, 32};
  int          halo = 4;
  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);

  // Tile (0, 0): inner x in [0, 32], inner y in [0, 32]
  TileRegion r00 = va.tile_region_from_tile_coords(0, 0);
  glm::vec4  inner00 = r00.inner_bbox();
  EXPECT_NEAR(inner00.x, 0.f, 1e-5f);
  EXPECT_NEAR(inner00.y, 32.f, 1e-5f);
  EXPECT_NEAR(inner00.z, 0.f, 1e-5f);
  EXPECT_NEAR(inner00.w, 32.f, 1e-5f);

  // Tile (1, 1): inner x in [32, 64], inner y in [32, 64]
  TileRegion r11 = va.tile_region_from_tile_coords(1, 1);
  glm::vec4  inner11 = r11.inner_bbox();
  EXPECT_NEAR(inner11.x, 32.f, 1e-5f);
  EXPECT_NEAR(inner11.y, 64.f, 1e-5f);
  EXPECT_NEAR(inner11.z, 32.f, 1e-5f);
  EXPECT_NEAR(inner11.w, 64.f, 1e-5f);
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
                    StorageMode::VA_DISK_LRU_STREAM,
                    StorageMode::VA_DISK_LRU_FULL})
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

  // 2. DISK_LRU_STREAM holds only actively pinned tiles (0 unpinned lingering
  // after trim)
  {
    VirtualArray va(shape, tile_shape, 2, StorageMode::VA_DISK_LRU_STREAM);
    va.fill(2.0f, cm);
    // In sequential stream mode, after fill finishes each tile is unpinned and
    // evicted
    EXPECT_EQ(va.live_tile_count(), 0);

    va.trim_storage();
    EXPECT_EQ(va.live_tile_count(), 0);
    EXPECT_EQ(va.live_memory_bytes(), 0);
  }

  // 3. DISK_LRU_FULL holds tiles in RAM until trim() is called
  {
    VirtualArray va(shape, tile_shape, 2, StorageMode::VA_DISK_LRU_FULL);
    va.fill(3.0f, cm);
    EXPECT_EQ(va.live_tile_count(), 16);

    va.trim_storage();
    EXPECT_EQ(va.live_tile_count(), 0);
    EXPECT_EQ(va.live_memory_bytes(), 0);
  }
}

TEST(VirtualArrayTest, SyncOverlapBuffersOperations)
{
  glm::ivec2 shape{32, 32};
  glm::ivec2 tile_shape{16, 16};
  int        halo = 4;

  auto setup_va = [&]()
  {
    VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);
    ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};
    for_each_tile(
        va,
        [](Array &tile, const TileRegion &region)
        {
          float base_val = float(region.key.tx + region.key.ty * 10);
          tile = base_val;
        },
        cm);
    return va;
  };

  // Test Average / Mean
  for (auto op : {SyncOperation::Average, SyncOperation::Mean})
  {
    VirtualArray va = setup_va();
    va.sync_overlap_buffers(op);

    TileRegion r0 = va.tile_region_from_tile_coords(0, 0);
    TileRegion r1 = va.tile_region_from_tile_coords(1, 0);
    Array     &t0 = va.storage->get_tile(r0);
    Array     &t1 = va.storage->get_tile(r1);

    int pbuf = t0.shape.x - 2 * halo;
    EXPECT_FLOAT_EQ(t0(pbuf, 0), 0.5f);
    EXPECT_FLOAT_EQ(t1(0, 0), 0.5f);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }

  // Test CopyFirst
  {
    VirtualArray va = setup_va();
    va.sync_overlap_buffers(SyncOperation::CopyFirst);

    TileRegion r0 = va.tile_region_from_tile_coords(0, 0);
    TileRegion r1 = va.tile_region_from_tile_coords(1, 0);
    Array     &t0 = va.storage->get_tile(r0);
    Array     &t1 = va.storage->get_tile(r1);

    int pbuf = t0.shape.x - 2 * halo;
    EXPECT_FLOAT_EQ(t0(pbuf, 0), 0.f);
    EXPECT_FLOAT_EQ(t1(0, 0), 0.f);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }

  // Test CopySecond
  {
    VirtualArray va = setup_va();
    va.sync_overlap_buffers(SyncOperation::CopySecond);

    TileRegion r0 = va.tile_region_from_tile_coords(0, 0);
    TileRegion r1 = va.tile_region_from_tile_coords(1, 0);
    Array     &t0 = va.storage->get_tile(r0);
    Array     &t1 = va.storage->get_tile(r1);

    int pbuf = t0.shape.x - 2 * halo;
    EXPECT_FLOAT_EQ(t0(pbuf, 0), 1.f);
    EXPECT_FLOAT_EQ(t1(0, 0), 1.f);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }

  // Test Min
  {
    VirtualArray va = setup_va();
    va.sync_overlap_buffers(SyncOperation::Min);

    TileRegion r0 = va.tile_region_from_tile_coords(0, 0);
    TileRegion r1 = va.tile_region_from_tile_coords(1, 0);
    Array     &t0 = va.storage->get_tile(r0);
    Array     &t1 = va.storage->get_tile(r1);

    int pbuf = t0.shape.x - 2 * halo;
    EXPECT_FLOAT_EQ(t0(pbuf, 0), 0.f);
    EXPECT_FLOAT_EQ(t1(0, 0), 0.f);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }

  // Test Max
  {
    VirtualArray va = setup_va();
    va.sync_overlap_buffers(SyncOperation::Max);

    TileRegion r0 = va.tile_region_from_tile_coords(0, 0);
    TileRegion r1 = va.tile_region_from_tile_coords(1, 0);
    Array     &t0 = va.storage->get_tile(r0);
    Array     &t1 = va.storage->get_tile(r1);

    int pbuf = t0.shape.x - 2 * halo;
    EXPECT_FLOAT_EQ(t0(pbuf, 0), 1.f);
    EXPECT_FLOAT_EQ(t1(0, 0), 1.f);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }

  // Test SmoothBlend
  {
    VirtualArray va = setup_va();
    va.sync_overlap_buffers(SyncOperation::SmoothBlend);

    TileRegion r0 = va.tile_region_from_tile_coords(0, 0);
    TileRegion r1 = va.tile_region_from_tile_coords(1, 0);
    Array     &t0 = va.storage->get_tile(r0);
    Array     &t1 = va.storage->get_tile(r1);

    int pbuf = t0.shape.x - 2 * halo;
    EXPECT_FLOAT_EQ(t0(pbuf, 0), t1(0, 0));
    EXPECT_FLOAT_EQ(t0(pbuf, 0), 0.f);

    int pbuf_end = t0.shape.x - halo - 1;
    EXPECT_FLOAT_EQ(t0(pbuf_end, 0), t1(halo - 1, 0));
    EXPECT_FLOAT_EQ(t0(pbuf_end, 0), 1.f);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }
}

TEST(VirtualArrayTest, SumPerTile)
{
  glm::ivec2   shape{64, 64};
  glm::ivec2   tile_shape{32, 32};
  int          halo = 4;
  VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};

  // Fill with 1.0f
  va.fill(1.0f, cm);

  std::vector<float> tile_sums;
  float              total_inner_sum = 0.f;

  sum_per_tile(va, tile_sums, &total_inner_sum, nullptr, 0.f, cm);

  // va has 2x2 tiles = 4 tiles total
  EXPECT_EQ(tile_sums.size(), 4u);

  // Total inner sum should equal 64 * 64 * 1.0f = 4096.0f
  EXPECT_NEAR(total_inner_sum, 64.f * 64.f, 1e-4f);

  // Each tile has shape 32 + 2 * halo (or boundary halo)
  for (float sum : tile_sums)
  {
    EXPECT_GT(sum, 0.f);
  }

  // Test with exclusion mask
  VirtualArray exclusion(shape, tile_shape, halo, StorageMode::VA_RAM);
  exclusion.fill(0.0f, cm);

  // Exclude top half
  Array excl_array(shape, 0.0f);
  for (int j = 32; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      excl_array(i, j) = 1.0f;
    }
  }
  exclusion.from_array(excl_array, cm);

  std::vector<float> tile_sums_excl;
  float              total_inner_sum_excl = 0.f;

  sum_per_tile(va, tile_sums_excl, &total_inner_sum_excl, &exclusion, 0.5f, cm);

  // Total inner sum should be half: 64 * 32 * 1.0f = 2048.0f
  EXPECT_NEAR(total_inner_sum_excl, 64.f * 32.f, 1e-4f);

  // Bottom row tiles have full inner contribution plus halo
  EXPECT_GT(tile_sums_excl[0], 0.f);
  EXPECT_GT(tile_sums_excl[1], 0.f);

  // Top row tiles (ty = 1) overlap with non-excluded row in halo (halo = 4),
  // but if we exclude row 32 - halo onwards, top tiles would be completely 0.
  // Check that top row tiles have strictly smaller sum than bottom row tiles.
  EXPECT_LT(tile_sums_excl[2], tile_sums_excl[0]);
  EXPECT_LT(tile_sums_excl[3], tile_sums_excl[1]);
}

TEST(VirtualArrayTest, MPICollectiveMinTiledVsNonTiled)
{
  glm::ivec2 shape{128, 128};
  glm::ivec2 tile_shape{32, 32};
  int        halo = 0;

  // 1. Reference monolithic Array
  Array full_array(shape, 100.0f);
  // Inject known minimum in a specific tile
  full_array(42, 75) = -17.25f;

  // Compute with default SingleArrayComm
  float     ref_min = full_array.min();
  float     ref_max = full_array.max();
  float     ref_sum = full_array.sum();
  float     ref_mean = full_array.mean();
  float     ref_ptp = full_array.ptp();
  glm::vec2 ref_range = full_array.range();
  float     ref_std = full_array.std();

  EXPECT_FLOAT_EQ(ref_min, -17.25f);
  EXPECT_FLOAT_EQ(ref_max, 100.0f);

  // 2. Test across compute mode: Sequential
  {
    VirtualArray va(shape, tile_shape, halo, StorageMode::VA_DISK_LRU_STREAM);
    ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL, .trim_storage = true};

    va.from_array(full_array, cm);

    float     va_min = std::numeric_limits<float>::infinity();
    float     va_max = -std::numeric_limits<float>::infinity();
    float     va_sum = 0.f;
    float     va_mean = 0.f;
    float     va_ptp = 0.f;
    glm::vec2 va_range{0.f, 0.f};

    for_each_tile(
        va,
        [&](const Array &tile, const TileRegion &, ICommunicator &comm)
        {
          va_min = tile.min(comm);
          va_max = tile.max(comm);
          va_sum = tile.sum(comm);
          va_mean = tile.mean(comm);
          va_ptp = tile.ptp(comm);
          va_range = tile.range(comm);
        },
        cm);

    EXPECT_FLOAT_EQ(va_min, ref_min);
    EXPECT_FLOAT_EQ(va_max, ref_max);
    EXPECT_NEAR(va_sum, ref_sum, 1e-3f);
    EXPECT_NEAR(va_mean, ref_mean, 1e-4f);
    EXPECT_FLOAT_EQ(va_ptp, ref_ptp);
    EXPECT_FLOAT_EQ(va_range.x, ref_range.x);
    EXPECT_FLOAT_EQ(va_range.y, ref_range.y);

    // Ensure bounded memory footprint under sequential execution
    EXPECT_LE(va.live_tile_count(), 2);
  }

  // 3. Test with SingleArrayComm (default standalone argument)
  {
    EXPECT_FLOAT_EQ(full_array.min(), ref_min);
    EXPECT_FLOAT_EQ(full_array.max(), ref_max);
    EXPECT_NEAR(full_array.sum(), ref_sum, 1e-3f);
    EXPECT_NEAR(full_array.mean(), ref_mean, 1e-4f);
    EXPECT_FLOAT_EQ(full_array.ptp(), ref_ptp);
    EXPECT_FLOAT_EQ(full_array.range().x, ref_range.x);
    EXPECT_FLOAT_EQ(full_array.range().y, ref_range.y);
    EXPECT_NEAR(full_array.std(), ref_std, 1e-4f);
  }

  // 4. Test unique_values() with communicator
  {
    Array label_array(shape, 1.0f);
    label_array(10, 10) = 2.0f;
    label_array(40, 50) = 5.0f;
    label_array(100, 110) = 9.0f;

    std::vector<float> ref_unique = label_array.unique_values();
    std::vector<float> expected{1.0f, 2.0f, 5.0f, 9.0f};
    EXPECT_EQ(ref_unique, expected);

    VirtualArray va_labels(shape,
                           tile_shape,
                           halo,
                           StorageMode::VA_DISK_LRU_STREAM);
    ComputeMode  cm_seq{.mode = ForEachMode::VA_SEQUENTIAL,
                        .trim_storage = true};
    va_labels.from_array(label_array, cm_seq);

    std::vector<float> va_unique;
    for_each_tile(
        va_labels,
        [&](const Array &tile, const TileRegion &, ICommunicator &comm)
        { va_unique = tile.unique_values(comm); },
        cm_seq);

    EXPECT_EQ(va_unique, ref_unique);
    EXPECT_LE(va_labels.live_tile_count(), 2);
  }

  // 5. Test sync_overlap_buffers() with communicator
  {
    int          va_halo = 4;
    VirtualArray va_sync(shape, tile_shape, va_halo, StorageMode::VA_RAM);
    ComputeMode  cm_seq{.mode = ForEachMode::VA_SEQUENTIAL};

    // Fill each tile with a value corresponding to its tile index
    for_each_tile(
        va_sync,
        [](Array &tile, const TileRegion &region, ICommunicator &comm)
        {
          float tile_val = static_cast<float>(region.key.tx +
                                              region.key.ty * 10);
          for (float &v : tile.vector)
            v = tile_val;

          // Call communicator sync_overlap_buffers directly
          comm.sync_overlap_buffers(tile, region, SyncOperation::CopyFirst);
        },
        cm_seq);

    // Verify west-to-east neighbor halo on tile (1, 0) received tile (0, 0)'s
    // value
    TileRegion   r10 = va_sync.tile_region_from_tile_coords(1, 0);
    const Array &tile10 = va_sync.storage->get_tile(r10);
    // Western halo is at x in [0, va_halo)
    EXPECT_FLOAT_EQ(tile10(0, tile10.shape.y / 2),
                    0.0f); // value from tile (0, 0)
    va_sync.storage->release_tile(r10);
  }

  // 6. Test comparison between VirtualArray post-sync and Communicator sync
  {
    const glm::ivec2 test_shape{128, 128};
    const glm::ivec2 test_tile{32, 32};
    const int        test_halo = 8;
    const int        smooth_ir = 16;

    Array base_arr = noise_fbm(NoiseType::PERLIN, test_shape, {4.f, 4.f}, 42);

    // Reference: Monolithic smooth
    Array ref_arr = base_arr;
    smooth_cpulse(ref_arr, smooth_ir);

    // Mode A: Post-sync via VirtualArray::sync_overlap_buffers()
    VirtualArray va_post(test_shape, test_tile, test_halo, StorageMode::VA_RAM);
    ComputeMode  cm_dist{.mode = ForEachMode::VA_DISTRIBUTED};
    va_post.from_array(base_arr, cm_dist);

    for_each_tile(
        va_post,
        [&](Array &tile, const TileRegion &)
        { smooth_cpulse(tile, smooth_ir); },
        cm_dist);
    va_post.sync_overlap_buffers();
    Array res_post = va_post.to_array(cm_dist);

    // Mode B: Communicator sync inside for_each_tile()
    VirtualArray va_comm(test_shape, test_tile, test_halo, StorageMode::VA_RAM);
    va_comm.from_array(base_arr, cm_dist);

    for_each_tile(
        va_comm,
        [&](Array &tile, const TileRegion &region, ICommunicator &comm)
        { gpu::smooth_cpulse(tile, smooth_ir, nullptr, region, comm); },
        cm_dist);
    Array res_comm = va_comm.to_array(cm_dist);

    // Both should significantly match each other across boundary regions
    Array diff_post_comm = hmap::abs(res_post - res_comm);
    EXPECT_LT(diff_post_comm.mean(), 0.01f);

    // Both should be close to monolithic reference
    Array diff_comm_ref = hmap::abs(res_comm - ref_arr);
    EXPECT_LT(diff_comm_ref.mean(), 0.05f);
  }

  // 7. Test remap and normalize_and_clip with communicator
  {
    Array full_test(shape);
    for (int j = 0; j < shape.y; ++j)
    {
      for (int i = 0; i < shape.x; ++i)
      {
        full_test(i, j) = static_cast<float>(i + 2 * j);
      }
    }

    // Monolithic reference results
    Array ref_remap = full_test;
    remap(ref_remap, 0.f, 1.f);

    Array ref_norm = full_test;
    normalize_and_clip(ref_norm, 2.0f);

    // Distributed test via VirtualArray
    VirtualArray va_remap(shape, tile_shape, 0, StorageMode::VA_RAM);
    ComputeMode  cm_dist{.mode = ForEachMode::VA_DISTRIBUTED};
    va_remap.from_array(full_test, cm_dist);

    for_each_tile(
        va_remap,
        [](Array &tile, const TileRegion &, ICommunicator &comm)
        { remap(tile, 0.f, 1.f, comm); },
        cm_dist);

    Array res_remap = va_remap.to_array(cm_dist);
    Array diff_remap = hmap::abs(res_remap - ref_remap);
    EXPECT_LT(diff_remap.max(), 1e-5f);

    VirtualArray va_norm(shape, tile_shape, 0, StorageMode::VA_RAM);
    va_norm.from_array(full_test, cm_dist);

    for_each_tile(
        va_norm,
        [](Array &tile, const TileRegion &, ICommunicator &comm)
        { normalize_and_clip(tile, 2.0f, true, comm); },
        cm_dist);

    Array res_norm = va_norm.to_array(cm_dist);
    Array diff_norm = hmap::abs(res_norm - ref_norm);
    EXPECT_LT(diff_norm.max(), 1e-4f);
  }
}

TEST(VirtualArrayTest, StorageClonePersistenceAndDeepCopy)
{
  glm::ivec2 shape{32, 32};
  glm::ivec2 tile_shape{16, 16};
  int        halo = 2;

  for (auto mode :
       {StorageMode::VA_DISK_LRU_STREAM, StorageMode::VA_DISK_LRU_FULL})
  {
    VirtualArray va(shape, tile_shape, halo, mode);
    ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};
    va.fill(42.5f, cm);

    // Deep clone
    auto cloned = va.clone(cm, true);
    EXPECT_FLOAT_EQ(cloned->min(cm), 42.5f);
    EXPECT_FLOAT_EQ(cloned->max(cm), 42.5f);

    // Modify original
    va.fill(10.0f, cm);
    EXPECT_FLOAT_EQ(va.min(cm), 10.0f);

    // Cloned must retain its own copy
    EXPECT_FLOAT_EQ(cloned->min(cm), 42.5f);
  }
}

TEST(VirtualArrayTest, SchedulerCooperativeExceptionPropagation)
{
  glm::ivec2   shape{32, 32};
  glm::ivec2   tile_shape{16, 16};
  VirtualArray va(shape, tile_shape, 2, StorageMode::VA_DISK_LRU_STREAM);
  ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};

  // One tile throws an exception inside collective operation
  auto faulty_lambda =
      [](Array &, const TileRegion &region, ICommunicator &comm)
  {
    if (region.key.tx == 1 && region.key.ty == 1)
    {
      throw std::runtime_error("intentional error on tile 1,1");
    }
    comm.barrier();
  };

  EXPECT_THROW(for_each_tile(va, faulty_lambda, cm), std::runtime_error);
}
