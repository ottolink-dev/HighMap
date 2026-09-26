/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <cmath>
#include <iostream>

#include "highmap.hpp"

using namespace hmap;

int main(void)
{
  // --- Create and synthesize a VirtualArray heightmap
  const glm::ivec2 total_shape = {1024, 1024};
  const glm::ivec2 tile_shape = {256, 256};
  const int        halo = 4;

  VirtualArray va(total_shape,
                  {-1.f, 1.f, -1.f, 1.f},
                  tile_shape,
                  halo,
                  StorageMode::VA_RAM);

  ComputeMode cm;
  cm.mode = ForEachMode::VA_DISTRIBUTED;

  // populate with a synthetic terrain pattern
  Array base_noise(total_shape);

  for (int j = 0; j < total_shape.y; ++j)
    for (int i = 0; i < total_shape.x; ++i)
    {
      base_noise(i, j) = 100.f + 50.f * std::sin(float(i) * 0.02f) *
                                     std::cos(float(j) * 0.02f);
    }

  va.from_array(base_noise, cm);

  // --- 1. OpenEXR Export (OpenEXRWriter)
  // Truly tiled and distributed streaming: each tile is written
  // directly to disk via Imf::TiledOutputFile with constant O(Tile)
  // memory footprint.
  {
    ImageWriterConfig cfg;
    cfg.data_type = ImageDataType::FLOAT32;
    cfg.image_shape = va.shape;
    cfg.tile_shape = va.tile_shape;
    cfg.channels = 1;
    cfg.enable_remapping = true;
    cfg.remap_min = 0.f;
    cfg.remap_max = 200.f;

    std::cout << "Exporting to OpenEXR (tiled streaming)..." << std::endl;
    export_image(va, "heightmap_tiled.exr", cfg, cm);
  }

  // --- 2. BigTIFF Export (BigTiffWriter)
  // Native tiled BigTIFF writing via libtiff (TIFFWriteTile).
  // Supports large rasters (> 4GB) with constant O(Tile) memory
  // footprint.
  {
    ImageWriterConfig cfg;
    cfg.data_type = ImageDataType::UINT16;
    cfg.image_shape = va.shape;
    cfg.tile_shape = va.tile_shape;
    cfg.channels = 1;
    cfg.compression = "deflate"; // or "lzw", "none"
    cfg.enable_remapping = true;
    cfg.remap_min = 0.f;
    cfg.remap_max = 200.f;

    std::cout << "Exporting to BigTIFF (tiled streaming)..." << std::endl;
    export_image(va, "heightmap_tiled.tif", cfg, cm);
  }

  // --- 4. OpenCV Export (OpenCVWriter: PNG, JPG, BMP)

  // NOTE & CAUTION: Unlike OpenEXR, BigTIFF, and Zarr, standard image
  // formats like PNG and JPEG do NOT natively support random-access
  // tiled disk writing.  OpenCVWriter acts as an in-memory
  // accumulator: tiles are assembled into an internal cv::Mat buffer
  // in RAM, and encoded via cv::imwrite upon close().  This is NOT
  // streaming or distributed/memory-bounded for huge rasters
  // exceeding system RAM.
  {
    ImageWriterConfig cfg;
    cfg.data_type = ImageDataType::UINT16; // 16-bit PNG
    cfg.image_shape = va.shape;
    cfg.tile_shape = va.tile_shape;
    cfg.channels = 1;
    cfg.enable_remapping = true;
    cfg.remap_min = 0.f;
    cfg.remap_max = 200.f;

    std::cout << "Exporting to 16-bit PNG via OpenCV (in-memory buffered)..."
              << std::endl;
    export_image(va, "heightmap_16bit.png", cfg, cm);
  }

  std::cout << "All exports completed successfully!" << std::endl;
  return 0;
}
