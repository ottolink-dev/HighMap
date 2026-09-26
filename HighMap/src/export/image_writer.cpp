/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <filesystem>
#include <memory>
#include <string>

#include "highmap/export/bigtiff_writer.hpp"
#include "highmap/export/image_writer.hpp"
#include "highmap/export/opencv_writer.hpp"
#include "highmap/export/openexr_writer.hpp"
#include "highmap/internal/string_utils.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

// =====================================
// ImageWriter default implementations
// =====================================

bool ImageWriter::write_chunk(int                       offset_x,
                              int                       offset_y,
                              const std::vector<Array> &channels)
{
  if (channels.empty()) return false;
  return this->write_chunk(offset_x, offset_y, channels[0]);
}

bool ImageWriter::write_tile(const TileRegion &region, const Array &data)
{
  int core_x = region.key.tx * this->config().tile_shape.x;
  int core_y = region.key.ty * this->config().tile_shape.y;
  int core_w = region.shape.x - region.halo.x - region.halo.y;
  int core_h = region.shape.y - region.halo.z - region.halo.w;

  if (core_w <= 0 || core_h <= 0) return false;

  if (region.halo.x > 0 || region.halo.y > 0 || region.halo.z > 0 ||
      region.halo.w > 0)
  {
    Array core_tile = data.extract_slice(region.halo.x,
                                         region.halo.x + core_w,
                                         region.halo.z,
                                         region.halo.z + core_h);
    return this->write_chunk(core_x, core_y, core_tile);
  }

  return this->write_chunk(core_x, core_y, data);
}

std::unique_ptr<ImageWriter> ImageWriter::create(
    const std::string       &filepath,
    const ImageWriterConfig &config)
{
  std::string ext = std::filesystem::path(filepath).extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

  std::unique_ptr<ImageWriter> writer;

  if (ext == ".exr")
  {
    writer = std::make_unique<OpenEXRWriter>(filepath, config);
  }
  else if (ext == ".tif" || ext == ".tiff" || ext == ".btif" ||
           ext == ".bigtiff")
  {
    writer = std::make_unique<BigTiffWriter>(filepath, config);
  }
  else
  {
    writer = std::make_unique<OpenCVWriter>(filepath, config);
  }

  return writer;
}

} // namespace hmap
