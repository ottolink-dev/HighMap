#include <filesystem>
#include <fstream>

#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/export.hpp"
#include "highmap/primitives.hpp"
#include "highmap/virtual_array/virtual_array.hpp"
#include "highmap/virtual_array/virtual_texture.hpp"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

using namespace hmap;

TEST(ImageWriterTest, FactoryCreation)
{
  ImageWriterConfig config;
  config.image_shape = {128, 128};
  config.tile_shape = {64, 64};

  auto exr_writer = ImageWriter::create("test.exr", config);
  EXPECT_NE(dynamic_cast<OpenEXRWriter *>(exr_writer.get()), nullptr);

  auto tiff_writer = ImageWriter::create("test.tif", config);
  EXPECT_NE(dynamic_cast<BigTiffWriter *>(tiff_writer.get()), nullptr);

  auto png_writer = ImageWriter::create("test.png", config);
  EXPECT_NE(dynamic_cast<OpenCVWriter *>(png_writer.get()), nullptr);

  exr_writer->close();
  tiff_writer->close();
  png_writer->close();

  if (std::filesystem::exists("test.exr")) std::filesystem::remove("test.exr");
  if (std::filesystem::exists("test.tif")) std::filesystem::remove("test.tif");
  if (std::filesystem::exists("test.png")) std::filesystem::remove("test.png");
}

TEST(ImageWriterTest, OpenEXRTiledRoundTrip)
{
  const std::string fname = "tmp_tiled_test.exr";
  if (std::filesystem::exists(fname)) std::filesystem::remove(fname);

  glm::ivec2 img_shape = {128, 128};
  glm::ivec2 tile_shape = {64, 64};

  Array full_array(img_shape);
  for (int j = 0; j < img_shape.y; ++j)
    for (int i = 0; i < img_shape.x; ++i)
      full_array(i, j) = static_cast<float>(i + j * img_shape.x);

  ImageWriterConfig config;
  config.image_shape = img_shape;
  config.tile_shape = tile_shape;
  config.data_type = ImageDataType::FLOAT32;

  OpenEXRWriter writer;
  ASSERT_TRUE(writer.open(fname, config));

  for (int ty = 0; ty < 2; ++ty)
  {
    for (int tx = 0; tx < 2; ++tx)
    {
      int   ox = tx * tile_shape.x;
      int   oy = ty * tile_shape.y;
      Array tile = full_array.extract_slice(ox,
                                            ox + tile_shape.x,
                                            oy,
                                            oy + tile_shape.y);
      EXPECT_TRUE(writer.write_chunk(ox, oy, tile));
    }
  }
  writer.close();

  ASSERT_TRUE(std::filesystem::exists(fname));

  Array loaded(fname, false);
  EXPECT_EQ(loaded.shape.x, img_shape.x);
  EXPECT_EQ(loaded.shape.y, img_shape.y);
  EXPECT_TRUE(assert_almost_equal(full_array, loaded, 1e-4f));

  std::filesystem::remove(fname);
}

TEST(ImageWriterTest, BigTiffTiledWrite)
{
  const std::string fname = "tmp_tiled_test.tif";
  if (std::filesystem::exists(fname)) std::filesystem::remove(fname);

  glm::ivec2 img_shape = {128, 128};
  glm::ivec2 tile_shape = {64, 64};

  Array full_array = white(img_shape, 0.f, 1.f, 42);

  ImageWriterConfig config;
  config.image_shape = img_shape;
  config.tile_shape = tile_shape;
  config.data_type = ImageDataType::FLOAT32;

  BigTiffWriter writer;
  ASSERT_TRUE(writer.open(fname, config));

  for (int ty = 0; ty < 2; ++ty)
  {
    for (int tx = 0; tx < 2; ++tx)
    {
      int   ox = tx * tile_shape.x;
      int   oy = ty * tile_shape.y;
      Array tile = full_array.extract_slice(ox,
                                            ox + tile_shape.x,
                                            oy,
                                            oy + tile_shape.y);
      EXPECT_TRUE(writer.write_chunk(ox, oy, tile));
    }
  }
  writer.close();

  ASSERT_TRUE(std::filesystem::exists(fname));

  Array loaded = read_to_array(fname, false, false);
  EXPECT_EQ(loaded.shape.x, img_shape.x);
  EXPECT_EQ(loaded.shape.y, img_shape.y);

  std::filesystem::remove(fname);
}

TEST(ImageWriterTest, VirtualArrayExportImage)
{
  const std::string fname = "tmp_va_export.exr";
  if (std::filesystem::exists(fname)) std::filesystem::remove(fname);

  glm::ivec2 shape = {128, 128};
  glm::ivec2 tile_shape = {64, 64};
  int        halo = 4;

  Array src(shape, 42.5f);

  VirtualArray va(shape,
                  {-1.f, 1.f, -1.f, 1.f},
                  tile_shape,
                  halo,
                  StorageMode::VA_RAM);
  ComputeMode  cm;
  cm.mode = ForEachMode::VA_DISTRIBUTED;

  va.fill(42.5f, cm);

  ImageWriterConfig cfg;
  cfg.data_type = ImageDataType::FLOAT32;

  EXPECT_TRUE(export_image(va, fname, cfg, cm));
  ASSERT_TRUE(std::filesystem::exists(fname));

  Array loaded(fname, false);
  EXPECT_EQ(loaded.shape.x, shape.x);
  EXPECT_EQ(loaded.shape.y, shape.y);
  EXPECT_NEAR(loaded(0, 0), 42.5f, 1e-4f);
  EXPECT_NEAR(loaded(100, 100), 42.5f, 1e-4f);

  std::filesystem::remove(fname);
}

TEST(ImageWriterTest, ExportImageVirtualTexture)
{
  std::string fname = "test_export_vt.exr";
  if (std::filesystem::exists(fname)) std::filesystem::remove(fname);

  glm::ivec2 shape = {128, 128};
  glm::ivec2 tile_shape = {64, 64};
  int        halo = 4;
  int        channels = 3;

  VirtualTexture vt(shape,
                    {-1.f, 1.f, -1.f, 1.f},
                    tile_shape,
                    halo,
                    channels,
                    StorageMode::VA_RAM);
  ComputeMode    cm;
  cm.mode = ForEachMode::VA_DISTRIBUTED;

  vt.fill(10.0f, cm);

  ImageWriterConfig cfg;
  cfg.data_type = ImageDataType::FLOAT32;

  EXPECT_TRUE(export_image(vt, fname, cfg, cm));
  ASSERT_TRUE(std::filesystem::exists(fname));

  std::filesystem::remove(fname);
}
