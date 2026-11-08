#include "highmap/dbg/assert.hpp"
#include "highmap/virtual_array/virtual_texture.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(VirtualTextureTest, ConstructionAndChannels)
{
  glm::ivec2 shape{48, 48};
  glm::ivec2 tile_shape{16, 16};
  int        halo = 2;

  VirtualTexture tex(shape, tile_shape, halo, 4, StorageMode::VA_RAM);

  EXPECT_EQ(tex.shape, shape);
  EXPECT_EQ(tex.tile_shape, tile_shape);
  EXPECT_EQ(tex.halo, halo);
  EXPECT_EQ(tex.channels(), 4);
  EXPECT_EQ(tex.channels_ptr().size(), 4);

  ComputeMode cm{.mode = ForEachMode::VA_SEQUENTIAL};
  tex.fill(1.0f, cm);

  for (int c = 0; c < 4; ++c)
  {
    EXPECT_FLOAT_EQ(tex.channel(c).mean(cm), 1.0f);
  }

  tex.fill(3, 0.5f, cm);
  EXPECT_FLOAT_EQ(tex.channel(3).mean(cm), 0.5f);
}

TEST(VirtualTextureTest, TextureRoundTrip)
{
  glm::ivec2     shape{32, 32};
  VirtualTexture tex(shape, {16, 16}, 2, 3, StorageMode::VA_RAM);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};

  tex.channel(0).fill(0.1f, cm);
  tex.channel(1).fill(0.5f, cm);
  tex.channel(2).fill(0.9f, cm);

  Texture t = tex.to_texture(shape, cm);
  EXPECT_EQ(t.shape, shape);
  EXPECT_EQ(t.num_channels(), 3);
  EXPECT_NEAR(t[0].mean(), 0.1f, 1e-5f);
  EXPECT_NEAR(t[1].mean(), 0.5f, 1e-5f);
  EXPECT_NEAR(t[2].mean(), 0.9f, 1e-5f);
}

TEST(VirtualTextureTest, ChannelConversion)
{
  glm::ivec2     shape{32, 32};
  VirtualTexture tex3(shape, {16, 16}, 2, 3, StorageMode::VA_RAM);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};
  tex3.fill(0.25f, cm);

  VirtualTexture tex4 = convert_texture_channels(tex3, 4, 1.0f, cm);
  EXPECT_EQ(tex4.channels(), 4);
  EXPECT_FLOAT_EQ(tex4.channel(0).mean(cm), 0.25f);
  EXPECT_FLOAT_EQ(tex4.channel(1).mean(cm), 0.25f);
  EXPECT_FLOAT_EQ(tex4.channel(2).mean(cm), 0.25f);
  EXPECT_FLOAT_EQ(tex4.channel(3).mean(cm), 1.0f);
}

TEST(VirtualTextureTest, Luminance)
{
  glm::ivec2     shape{32, 32};
  VirtualTexture tex(shape, {16, 16}, 2, 3, StorageMode::VA_RAM);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};

  tex.channel(0).fill(1.0f, cm); // R = 1.0
  tex.channel(1).fill(0.0f, cm); // G = 0.0
  tex.channel(2).fill(0.0f, cm); // B = 0.0

  VirtualArray lum(shape, {16, 16}, 2, StorageMode::VA_RAM);
  luminance(lum, tex, cm);

  // Expected perceptual luminance = 0.299 * 1.0 = 0.299
  EXPECT_NEAR(lum.mean(cm), 0.299f, 1e-5f);
}

TEST(VirtualTextureTest, MixFourChannels)
{
  glm::ivec2     shape{32, 32};
  VirtualTexture t1(shape, {16, 16}, 2, 4, StorageMode::VA_RAM);
  VirtualTexture t2(shape, {16, 16}, 2, 4, StorageMode::VA_RAM);
  VirtualTexture out(shape, {16, 16}, 2, 4, StorageMode::VA_RAM);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};

  t1.fill(0.0f, cm);
  t1.fill(3, 1.0f, cm); // Alpha = 1.0

  t2.fill(1.0f, cm);
  t2.fill(3, 1.0f, cm); // Alpha = 1.0

  mix(out, t1, t2, cm, MixMethod::MM_LINEAR);

  EXPECT_NEAR(out.channel(0).mean(cm), 1.0f, 1e-5f);
  EXPECT_NEAR(out.channel(3).mean(cm), 1.0f, 1e-5f);
}

TEST(VirtualTextureTest, ColorizeThreeChannelsDoesNotCrash)
{
  glm::ivec2     shape{32, 32};
  VirtualArray   level(shape, {16, 16}, 2, StorageMode::VA_RAM);
  VirtualTexture tex3(shape, {16, 16}, 2, 3, StorageMode::VA_RAM);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};

  level.fill(0.5f, cm);

  // This previously crashed due to unconditional write to p_arrays[6] (alpha)
  EXPECT_NO_THROW({ colorize(tex3, level, cm, 0.f, 1.f, Cmap::VIRIDIS); });

  EXPECT_GT(tex3.channel(0).mean(cm), 0.f);
}

TEST(VirtualTextureTest, MixNormalMapThreeChannels)
{
  glm::ivec2     shape{32, 32};
  VirtualTexture base(shape, {16, 16}, 2, 3, StorageMode::VA_RAM);
  VirtualTexture detail(shape, {16, 16}, 2, 3, StorageMode::VA_RAM);
  VirtualTexture out(shape, {16, 16}, 2, 3, StorageMode::VA_RAM);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};

  base.channel(0).fill(0.5f, cm);
  base.channel(1).fill(0.5f, cm);
  base.channel(2).fill(1.0f, cm);

  detail.channel(0).fill(0.5f, cm);
  detail.channel(1).fill(0.5f, cm);
  detail.channel(2).fill(1.0f, cm);

  EXPECT_NO_THROW({
    mix_normal_map(out,
                   base,
                   detail,
                   cm,
                   1.0f,
                   NormalMapBlendingMethod::NMAP_LINEAR);
  });

  // Flat normal mapped back to RGB [0.5, 0.5, 1.0]
  EXPECT_NEAR(out.channel(0).mean(cm), 0.5f, 1e-3f);
  EXPECT_NEAR(out.channel(1).mean(cm), 0.5f, 1e-3f);
  EXPECT_NEAR(out.channel(2).mean(cm), 1.0f, 1e-3f);
}

TEST(VirtualTextureTest, MemoryFootprintAndTrim)
{
  glm::ivec2     shape{32, 32};
  VirtualTexture tex(shape, {16, 16}, 2, 4, StorageMode::VA_DISK_LRU_FULL);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};

  tex.fill(1.0f, cm);
  EXPECT_GT(tex.live_memory_bytes(), 0);

  tex.trim_storage();
  EXPECT_EQ(tex.live_memory_bytes(), 0);
}

TEST(VirtualTextureTest, SyncOverlapBuffers)
{
  glm::ivec2     shape{32, 16};
  glm::ivec2     tile_shape{16, 16};
  int            halo = 2;
  VirtualTexture tex(shape, tile_shape, halo, 2, StorageMode::VA_RAM);
  ComputeMode    cm{.mode = ForEachMode::VA_SEQUENTIAL};

  // Fill channel 0 tile 0 with 0.0f, tile 1 with 1.0f
  // Fill channel 1 tile 0 with 2.0f, tile 1 with 4.0f
  for_each_tile(
      tex,
      [&](std::vector<Array *> &chs, const TileRegion &region)
      {
        if (region.key.tx == 0)
        {
          *chs[0] = 0.0f;
          *chs[1] = 2.0f;
        }
        else
        {
          *chs[0] = 1.0f;
          *chs[1] = 4.0f;
        }
      },
      cm);

  tex.sync_overlap_buffers(SyncOperation::Average);

  for (int c = 0; c < tex.channels(); ++c)
  {
    auto      &va = tex.channel(c);
    TileRegion r0 = va.tile_region_from_tile_coords(0, 0);
    TileRegion r1 = va.tile_region_from_tile_coords(1, 0);
    Array     &t0 = va.storage->get_tile(r0);
    Array     &t1 = va.storage->get_tile(r1);

    int   pbuf = t0.shape.x - 2 * halo;
    float expected = (c == 0) ? 0.5f : 3.0f;

    EXPECT_FLOAT_EQ(t0(pbuf, 0), expected);
    EXPECT_FLOAT_EQ(t1(0, 0), expected);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }
}
