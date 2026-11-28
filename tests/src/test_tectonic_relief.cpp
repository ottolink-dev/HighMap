/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#include "highmap/primitives/geo.hpp"

#include <gtest/gtest.h>

TEST(TectonicRelief, ConvergenceVsDivergence)
{
  const glm::ivec2 shape = {64, 64};
  const float      roughness = 0.5f;
  const float      boundary_width = 0.15f;

  // two continental plates converging toward x = 0.5
  std::vector<hmap::TectonicPlate> conv_plates = {
      {{0.25f, 0.5f}, {0.8f, 0.0f}, true},
      {{0.75f, 0.5f}, {-0.8f, 0.0f}, true},
  };

  hmap::Array conv_ranges;
  hmap::Array conv_h = hmap::tectonic_relief(shape,
                                             conv_plates,
                                             roughness,
                                             boundary_width,
                                             1.0f,
                                             0.8f,
                                             0.15f,
                                             false,
                                             nullptr,
                                             nullptr,
                                             nullptr,
                                             nullptr,
                                             &conv_ranges);

  // two continental plates diverging from x = 0.5
  std::vector<hmap::TectonicPlate> div_plates = {
      {{0.25f, 0.5f}, {-0.8f, 0.0f}, true},
      {{0.75f, 0.5f}, {0.8f, 0.0f}, true},
  };

  hmap::Array div_ranges;
  hmap::Array div_h = hmap::tectonic_relief(shape,
                                            div_plates,
                                            roughness,
                                            boundary_width,
                                            1.0f,
                                            0.8f,
                                            0.15f,
                                            false,
                                            nullptr,
                                            nullptr,
                                            nullptr,
                                            nullptr,
                                            &div_ranges);

  int mid_x = shape.x / 2;
  int mid_y = shape.y / 2;

  EXPECT_GT(conv_h(mid_x, mid_y), div_h(mid_x, mid_y));
  EXPECT_GT(conv_ranges(mid_x, mid_y), 0.1f);
  EXPECT_FLOAT_EQ(div_ranges(mid_x, mid_y), 0.0f);
}

TEST(TectonicRelief, MultiPlateSeaLandBalance)
{
  const glm::ivec2 shape = {128, 128};
  hmap::Array      h = hmap::tectonic_relief(shape, 9, 1234, 0.42f);

  int below_sea = 0;
  int total = shape.x * shape.y;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if (h(i, j) < 0.0f) ++below_sea;

  float sea_ratio = static_cast<float>(below_sea) / static_cast<float>(total);

  EXPECT_GE(sea_ratio, 0.15f);
  EXPECT_LE(sea_ratio, 0.85f);
}

TEST(TectonicRelief, TileConsistency)
{
  const glm::ivec2 full_shape = {64, 64};
  const glm::ivec2 tile_shape = {32, 32};

  hmap::Array full_h = hmap::tectonic_relief(full_shape, 9, 99);

  // evaluate 2x2 tiles
  for (int ty = 0; ty < 2; ++ty)
  {
    for (int tx = 0; tx < 2; ++tx)
    {
      hmap::UvRect rect(tx * 0.5f, (tx + 1) * 0.5f, ty * 0.5f, (ty + 1) * 0.5f);
      hmap::Array  tile_h = hmap::tectonic_relief(tile_shape,
                                                 9,
                                                 99,
                                                 0.42f,
                                                 0.5f,
                                                 0.5f,
                                                 0.08f,
                                                 1.0f,
                                                 0.8f,
                                                 0.15f,
                                                 false,
                                                 nullptr,
                                                 nullptr,
                                                 nullptr,
                                                 nullptr,
                                                 nullptr,
                                                 nullptr,
                                                 nullptr,
                                                 nullptr,
                                                 rect);

      for (int j = 0; j < tile_shape.y; ++j)
      {
        for (int i = 0; i < tile_shape.x; ++i)
        {
          int gx = tx * tile_shape.x + i;
          int gy = ty * tile_shape.y + j;
          EXPECT_NEAR(tile_h(i, j), full_h(gx, gy), 1e-5f);
        }
      }
    }
  }
}
