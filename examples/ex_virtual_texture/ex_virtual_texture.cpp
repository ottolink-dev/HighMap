#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "highmap.hpp"

#define EXPECT(cond, msg)                                                      \
  do                                                                           \
  {                                                                            \
    if (!(cond))                                                               \
      std::cerr << "[FAIL] " << msg << "\n";                                   \
    else                                                                       \
      std::cout << "[ OK ] " << msg << "\n";                                   \
  } while (0)

int main()
{
  // ===========================================================================
  // Configuration
  // ===========================================================================

  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 shape{1000, 581};
  const glm::ivec2 tile_shape{220, 312};
  const size_t     halo = 32;

  const auto foreach_mode = hmap::ForEachMode::VA_DISTRIBUTED;
  const auto storage_mode = hmap::StorageMode::VA_DISK_LRU;

  const hmap::ComputeMode cm{
      .mode = foreach_mode,
      .trim_storage = false,
  };

  // ===========================================================================
  // Channel-based VirtualTexture
  // ===========================================================================

  hmap::VirtualTexture tex(shape, bbox, tile_shape, halo, 3, storage_mode);

  const auto fill_channels =
      [](std::vector<hmap::Array *> &tiles, const hmap::TileRegion &region)
  {
    auto &R = *tiles[0];
    auto &G = *tiles[1];
    auto &B = *tiles[2];

    for (int j = 0; j < region.shape.y; ++j)
      for (int i = 0; i < region.shape.x; ++i)
      {
        glm::vec2 pos = region.cell_center(i, j);
        R(i, j) = pos.x;
        G(i, j) = 0.f;
        B(i, j) = pos.y;
      }
  };

  hmap::for_each_tile(tex, fill_channels, cm);
  tex.to_png("out0.png", cm);

  // ===========================================================================
  // Tests: channel invariants
  // ===========================================================================

  {
    auto R = tex.channel(0).to_array(cm);
    auto G = tex.channel(1).to_array(cm);
    auto B = tex.channel(2).to_array(cm);

    EXPECT(R.shape == shape, "R channel shape");
    EXPECT(G.shape == shape, "G channel shape");
    EXPECT(B.shape == shape, "B channel shape");

    EXPECT(G.min() == 0.f, "G channel == 0");
  }

  // ===========================================================================
  // Per-pixel VirtualTexture
  // ===========================================================================

  hmap::VirtualTexture tex_px(shape, bbox, tile_shape, halo, 4, storage_mode);

  const auto fill_pixel =
      [](std::vector<float> &px, int i, int j, const hmap::TileRegion &region)
  {
    glm::vec2 pos = region.cell_center(i, j);
    px[0] = pos.x;
    px[1] = 0.f;
    px[2] = pos.y;
    px[3] = 1.f;
  };

  hmap::for_each_pixel(tex_px, fill_pixel, cm);
  tex_px.to_png("out1.png", cm);

  // ===========================================================================
  // Tests: per-pixel alpha
  // ===========================================================================

  {
    auto A = tex_px.channel(3).to_array(cm);

    EXPECT(A.min() == 1.f, "alpha min == 1");
    EXPECT(A.max() == 1.f, "alpha max == 1");
  }

  // ===========================================================================
  // Colorize test
  // ===========================================================================

  hmap::VirtualArray level(shape, bbox, tile_shape, halo, storage_mode);

  const auto gen_noise = [](hmap::Array &tile, const hmap::TileRegion &region)
  {
    tile = hmap::noise(hmap::NoiseType::PERLIN,
                       region.shape,
                       {2.f, 2.f},
                       0,
                       nullptr,
                       nullptr,
                       region.bbox);
  };

  hmap::for_each_tile(level, gen_noise, cm);
  level.remap(0.f, 1.f, cm);

  hmap::VirtualTexture tex_color(shape,
                                 bbox,
                                 tile_shape,
                                 halo,
                                 4,
                                 storage_mode);

  tex_color.fill(3, 1.f, cm);

  hmap::colorize(tex_color,
                 level,
                 cm,
                 level.min(cm),
                 level.max(cm),
                 hmap::Cmap::JET,
                 nullptr);

  tex_color.to_png("out2.png", cm);

  // ===========================================================================
  // Tests: colorize invariants
  // ===========================================================================

  {
    auto A = tex_color.channel(3).to_array(cm);
    auto R = tex_color.channel(0).to_array(cm);
    auto G = tex_color.channel(1).to_array(cm);
    auto B = tex_color.channel(2).to_array(cm);

    EXPECT(A.min() == 1.f, "colorize alpha preserved");
    EXPECT(R.max() <= 1.f, "R <= 1");
    EXPECT(G.max() <= 1.f, "G <= 1");
    EXPECT(B.max() <= 1.f, "B <= 1");
  }

  // ===========================================================================
  // Mix test
  // ===========================================================================

  hmap::VirtualTexture tex_mix(shape, bbox, tile_shape, halo, 4, storage_mode);

  hmap::mix(tex_mix, tex_px, tex_color, cm);
  tex_mix.to_png("out3.png", cm);

  // ===========================================================================
  // Tests: mix sanity
  // ===========================================================================

  {
    auto A = tex_mix.channel(3).to_array(cm);

    EXPECT(A.min() >= 0.f, "mix alpha >= 0");
    EXPECT(A.max() <= 1.f, "mix alpha <= 1");
  }

  // ===========================================================================
  // Small tile stress test
  // ===========================================================================

  {
    const glm::ivec2 small_tile{8, 8};
    const size_t     big_halo = 7;

    hmap::VirtualTexture tex_small(shape,
                                   bbox,
                                   small_tile,
                                   big_halo,
                                   3,
                                   storage_mode);

    hmap::for_each_tile(tex_small, fill_channels, cm);

    auto R = tex_small.channel(0).to_array(cm);
    auto B = tex_small.channel(2).to_array(cm);

    EXPECT(R.shape == shape, "small tile R shape");
    EXPECT(std::isfinite(R.min()), "small tile finite R");
    EXPECT(std::isfinite(B.max()), "small tile finite B");
  }

  return 0;
}
