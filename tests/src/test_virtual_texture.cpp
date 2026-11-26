#include "highmap/colormaps.hpp"
#include "highmap/internal/assert.hpp"
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
      [&](std::vector<Array *> &chs, const UvTile &region)
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
    auto  &va = tex.channel(c);
    UvTile r0 = va.tile_region_from_tile_coords(0, 0);
    UvTile r1 = va.tile_region_from_tile_coords(1, 0);
    Array &t0 = va.storage->get_tile(r0);
    Array &t1 = va.storage->get_tile(r1);

    int   pbuf = t0.shape.x - 2 * halo;
    float expected = (c == 0) ? 0.5f : 3.0f;

    EXPECT_FLOAT_EQ(t0(pbuf, 0), expected);
    EXPECT_FLOAT_EQ(t1(0, 0), expected);

    va.storage->release_tile(r0);
    va.storage->release_tile(r1);
  }
}

TEST(VirtualTextureTest, ColorizeTiledMatchesNonTiled)
{
  const glm::ivec2 shape{64, 64};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  Array level(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      level(i, j) = float(i + j * shape.x) / float(shape.x * shape.y - 1);

  Texture expected =
      colorize(level, 0.f, 1.f, Cmap::VIRIDIS, false, false, nullptr);

  ComputeMode cm{.mode = ForEachMode::VA_DISTRIBUTED};

  VirtualArray va_level(shape, tile_shape, halo, StorageMode::VA_RAM);
  va_level.from_array(level, cm);

  VirtualTexture va_tex(shape, tile_shape, halo, 3, StorageMode::VA_RAM);

  for_each_tile(
      {&va_level},
      {&va_tex.channel(0), &va_tex.channel(1), &va_tex.channel(2)},
      [](std::vector<const Array *> p_in,
         std::vector<Array *>       p_out,
         const UvTile &)
      {
        Texture tile_tex =
            colorize(*p_in[0], 0.f, 1.f, Cmap::VIRIDIS, false, false, nullptr);
        for (int c = 0; c < 3; ++c)
          *p_out[c] = tile_tex.channel(c);
      },
      cm);

  for (int c = 0; c < 3; ++c)
  {
    Array result = va_tex.channel(c).to_array(cm);
    EXPECT_TRUE(assert_almost_equal(result, expected.channel(c), 1e-5f));
  }
}

TEST(VirtualTextureTest, ColorizeBivariateTiledMatchesNonTiled)
{
  const glm::ivec2 shape{64, 64};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  Array a1(shape);
  Array a2(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      a1(i, j) = float(i) / float(shape.x - 1);
      a2(i, j) = float(j) / float(shape.y - 1);
    }

  const std::vector<float>     pos1{0.f, 1.f};
  const std::vector<float>     pos2{0.f, 1.f};
  const std::vector<glm::vec3> col1{{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}};
  const std::vector<glm::vec3> col2{{0.f, 0.f, 0.f}, {0.f, 1.f, 0.f}};

  Texture expected = colorize_bivariate(a1,
                                        a2,
                                        {0.f, 1.f},
                                        {0.f, 1.f},
                                        pos1,
                                        pos2,
                                        col1,
                                        col2,
                                        MixMethod::MM_SQRT_AVG,
                                        false,
                                        false,
                                        nullptr,
                                        nullptr);

  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};
  VirtualArray va_a1(shape, tile_shape, halo, StorageMode::VA_RAM);
  VirtualArray va_a2(shape, tile_shape, halo, StorageMode::VA_RAM);
  va_a1.from_array(a1, cm);
  va_a2.from_array(a2, cm);

  VirtualTexture va_tex(shape, tile_shape, halo, 3, StorageMode::VA_RAM);

  for_each_tile(
      {&va_a1, &va_a2},
      {&va_tex.channel(0), &va_tex.channel(1), &va_tex.channel(2)},
      [&](std::vector<const Array *> p_in,
          std::vector<Array *>       p_out,
          const UvTile &)
      {
        Texture tile_tex = colorize_bivariate(*p_in[0],
                                              *p_in[1],
                                              {0.f, 1.f},
                                              {0.f, 1.f},
                                              pos1,
                                              pos2,
                                              col1,
                                              col2,
                                              MixMethod::MM_SQRT_AVG,
                                              false,
                                              false,
                                              nullptr,
                                              nullptr);
        for (int c = 0; c < 3; ++c)
          *p_out[c] = tile_tex.channel(c);
      },
      cm);

  for (int c = 0; c < 3; ++c)
  {
    Array result = va_tex.channel(c).to_array(cm);
    EXPECT_TRUE(assert_almost_equal(result, expected.channel(c), 1e-5f));
  }
}

TEST(VirtualTextureTest, LuminanceTiledMatchesNonTiled)
{
  const glm::ivec2 shape{64, 64};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  Texture tex(shape, 3);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      tex.channel(0)(i, j) = float(i) / float(shape.x - 1);
      tex.channel(1)(i, j) = float(j) / float(shape.y - 1);
      tex.channel(2)(i, j) = 0.5f;
    }

  Array expected = luminance(tex);

  ComputeMode    cm{.mode = ForEachMode::VA_DISTRIBUTED};
  VirtualTexture va_tex(shape, tile_shape, halo, 3, StorageMode::VA_RAM);
  va_tex.from_arrays({&tex.channel(0), &tex.channel(1), &tex.channel(2)}, cm);

  VirtualArray va_lum(shape, tile_shape, halo, StorageMode::VA_RAM);

  for_each_tile(
      {&va_tex.channel(0), &va_tex.channel(1), &va_tex.channel(2)},
      {&va_lum},
      [](std::vector<const Array *> p_in,
         std::vector<Array *>       p_out,
         const UvTile &)
      {
        Texture tile_tex(p_in[0]->shape, 3);
        for (int c = 0; c < 3; ++c)
          tile_tex.channel(c) = *p_in[c];
        *p_out[0] = luminance(tile_tex);
      },
      cm);

  Array result = va_lum.to_array(cm);
  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-5f));
}

TEST(VirtualTextureTest, MixTiledMatchesNonTiled)
{
  const glm::ivec2 shape{64, 64};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  Texture t1(shape, 4);
  Texture t2(shape, 4);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float u = float(i) / float(shape.x - 1);
      float v = float(j) / float(shape.y - 1);
      t1.channel(0)(i, j) = u;
      t1.channel(1)(i, j) = 0.2f;
      t1.channel(2)(i, j) = 0.8f;
      t1.channel(3)(i, j) = 1.0f;
      t2.channel(0)(i, j) = 0.3f;
      t2.channel(1)(i, j) = v;
      t2.channel(2)(i, j) = 0.4f;
      t2.channel(3)(i, j) = 1.0f;
    }

  Texture expected = mix(t1, t2, MixMethod::MM_LINEAR);

  ComputeMode    cm{.mode = ForEachMode::VA_DISTRIBUTED};
  VirtualTexture va_t1(shape, tile_shape, halo, 4, StorageMode::VA_RAM);
  VirtualTexture va_t2(shape, tile_shape, halo, 4, StorageMode::VA_RAM);
  va_t1.from_arrays(
      {&t1.channel(0), &t1.channel(1), &t1.channel(2), &t1.channel(3)},
      cm);
  va_t2.from_arrays(
      {&t2.channel(0), &t2.channel(1), &t2.channel(2), &t2.channel(3)},
      cm);

  VirtualTexture va_out(shape, tile_shape, halo, 4, StorageMode::VA_RAM);

  for_each_tile(
      {&va_t1.channel(0),
       &va_t1.channel(1),
       &va_t1.channel(2),
       &va_t1.channel(3),
       &va_t2.channel(0),
       &va_t2.channel(1),
       &va_t2.channel(2),
       &va_t2.channel(3)},
      {&va_out.channel(0),
       &va_out.channel(1),
       &va_out.channel(2),
       &va_out.channel(3)},
      [](std::vector<const Array *> p_in,
         std::vector<Array *>       p_out,
         const UvTile &)
      {
        Texture tile_t1(p_in[0]->shape, 4);
        Texture tile_t2(p_in[0]->shape, 4);
        for (int c = 0; c < 4; ++c)
        {
          tile_t1.channel(c) = *p_in[c];
          tile_t2.channel(c) = *p_in[c + 4];
        }
        Texture tile_out = mix(tile_t1, tile_t2, MixMethod::MM_LINEAR);
        for (int c = 0; c < 4; ++c)
          *p_out[c] = tile_out.channel(c);
      },
      cm);

  for (int c = 0; c < 4; ++c)
  {
    Array result = va_out.channel(c).to_array(cm);
    EXPECT_TRUE(assert_almost_equal(result, expected.channel(c), 1e-5f));
  }
}

TEST(VirtualTextureTest, MixNormalMapTiledMatchesNonTiled)
{
  const glm::ivec2 shape{64, 64};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  Texture n_base(shape, 3);
  Texture n_detail(shape, 3);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      n_base.channel(0)(i, j) = 0.5f;
      n_base.channel(1)(i, j) = 0.5f;
      n_base.channel(2)(i, j) = 1.0f;
      n_detail.channel(0)(i, j) = 0.6f;
      n_detail.channel(1)(i, j) = 0.4f;
      n_detail.channel(2)(i, j) = 0.9f;
    }

  Texture expected = mix_normal_map(n_base,
                                    n_detail,
                                    0.5f,
                                    NormalMapBlendingMethod::NMAP_LINEAR);

  ComputeMode    cm{.mode = ForEachMode::VA_DISTRIBUTED};
  VirtualTexture va_base(shape, tile_shape, halo, 3, StorageMode::VA_RAM);
  VirtualTexture va_detail(shape, tile_shape, halo, 3, StorageMode::VA_RAM);
  va_base.from_arrays(
      {&n_base.channel(0), &n_base.channel(1), &n_base.channel(2)},
      cm);
  va_detail.from_arrays(
      {&n_detail.channel(0), &n_detail.channel(1), &n_detail.channel(2)},
      cm);

  VirtualTexture va_out(shape, tile_shape, halo, 3, StorageMode::VA_RAM);

  for_each_tile(
      {&va_base.channel(0),
       &va_base.channel(1),
       &va_base.channel(2),
       &va_detail.channel(0),
       &va_detail.channel(1),
       &va_detail.channel(2)},
      {&va_out.channel(0), &va_out.channel(1), &va_out.channel(2)},
      [](std::vector<const Array *> p_in,
         std::vector<Array *>       p_out,
         const UvTile &)
      {
        Texture tile_base(p_in[0]->shape, 3);
        Texture tile_detail(p_in[0]->shape, 3);
        for (int c = 0; c < 3; ++c)
        {
          tile_base.channel(c) = *p_in[c];
          tile_detail.channel(c) = *p_in[c + 3];
        }
        Texture tile_out = mix_normal_map(tile_base,
                                          tile_detail,
                                          0.5f,
                                          NormalMapBlendingMethod::NMAP_LINEAR);
        for (int c = 0; c < 3; ++c)
          *p_out[c] = tile_out.channel(c);
      },
      cm);

  for (int c = 0; c < 3; ++c)
  {
    Array result = va_out.channel(c).to_array(cm);
    EXPECT_TRUE(assert_almost_equal(result, expected.channel(c), 1e-5f));
  }
}
