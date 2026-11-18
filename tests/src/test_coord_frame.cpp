#include "highmap/array.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/interpolate/interpolate_terrain.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(CoordFrameTest, TransformAndIntersection)
{
  CoordFrame f1({0.f, 0.f}, {10.f, 10.f}, 0.f);
  CoordFrame f2({5.f, 5.f}, {10.f, 10.f}, 45.f);
  CoordFrame f3({100.f, 100.f}, {10.f, 10.f}, 0.f);

  EXPECT_TRUE(f1.intersects(f2));
  EXPECT_FALSE(f1.intersects(f3));

  glm::vec2 rel = f1.map_to_relative_coords(f2, 0.5f, 0.5f);
  glm::vec2 g1 = f1.map_to_global_coords(0.5f, 0.5f);
  glm::vec2 rel_expected = f2.map_to_relative_coords(g1.x, g1.y);

  EXPECT_NEAR(rel.x, rel_expected.x, 1e-5f);
  EXPECT_NEAR(rel.y, rel_expected.y, 1e-5f);
}

TEST(CoordFrameTest, ArrayIdentityPreservation)
{
  glm::ivec2 shape{32, 24};
  Array      src(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      src(i, j) = static_cast<float>(i * 3 + j * 7);

  CoordFrame frame({10.f, -5.f}, {50.f, 30.f}, 27.f);
  Array      dst(shape);

  interpolate_terrain(src, dst, frame, frame, /*bilinear=*/true);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      EXPECT_NEAR(dst(i, j), src(i, j), 1e-4f);
}

TEST(CoordFrameTest, VirtualArrayModesAndStreaming)
{
  glm::ivec2 shape{64, 64};
  glm::ivec2 tile_shape{16, 16};
  int        halo = 2;

  Array src_arr(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      src_arr(i, j) = static_cast<float>(i + 2 * j);

  CoordFrame t_src({0.f, 0.f}, {100.f, 100.f}, 0.f);
  CoordFrame t_dst({25.f, 25.f}, {50.f, 50.f}, 15.f);

  Array ref_dst(shape);
  interpolate_terrain(src_arr, ref_dst, t_src, t_dst, /*bilinear=*/true, -1.f);

  for (auto mode : {ForEachMode::VA_SEQUENTIAL,
                    ForEachMode::VA_DISTRIBUTED,
                    ForEachMode::VA_SINGLE_ARRAY})
  {
    StorageMode storage = (mode == ForEachMode::VA_SEQUENTIAL)
                              ? StorageMode::VA_DISK_LRU_STREAM
                              : StorageMode::VA_RAM;
    ComputeMode cm{.mode = mode};

    VirtualArray va_src(shape, tile_shape, halo, storage);
    va_src.from_array(src_arr, cm);

    VirtualArray va_dst(shape, tile_shape, halo, storage);
    interpolate_terrain(va_src,
                        va_dst,
                        t_src,
                        t_dst,
                        cm,
                        /*bilinear=*/true,
                        -1.f);

    if (storage == StorageMode::VA_DISK_LRU_STREAM)
    {
      EXPECT_EQ(va_src.live_tile_count(), 0u);
      EXPECT_EQ(va_dst.live_tile_count(), 0u);
    }

    Array dst_arr = va_dst.to_array(cm);
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
        EXPECT_NEAR(dst_arr(i, j), ref_dst(i, j), 1e-3f)
            << "Mismatch at (" << i << ", " << j << ") in mode "
            << to_string(mode);
  }
}

TEST(CoordFrameTest, FlattenHeightmapInPlace)
{
  glm::ivec2  shape{32, 32};
  glm::ivec2  tile_shape{16, 16};
  ComputeMode cm{.mode = ForEachMode::VA_DISTRIBUTED};

  VirtualArray va1(shape, tile_shape, 2, StorageMode::VA_RAM);
  VirtualArray va2(shape, tile_shape, 2, StorageMode::VA_RAM);
  va1.fill(1.0f, cm);
  va2.fill(5.0f, cm);

  CoordFrame t1({0.f, 0.f}, {10.f, 10.f}, 0.f);
  CoordFrame t2({2.f, 2.f}, {6.f, 6.f}, 0.f);

  flatten_terrain(va1, va2, t1, t2, cm);

  Array out = va1.to_array(cm);
  // Outside t2 remains 1.0, center of t2 blends toward 5.0
  EXPECT_FLOAT_EQ(out(0, 0), 1.0f);
  EXPECT_NEAR(out(16, 16), 5.0f, 0.1f);
}
