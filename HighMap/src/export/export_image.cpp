/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <memory>
#include <string>

#include "highmap/array.hpp"
#include "highmap/export/image.hpp"
#include "highmap/export/image_writer.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

bool export_image(const Array             &array,
                  const std::string       &fname,
                  const ImageWriterConfig &config)
{
  if (!validate_non_empty(array)) return false;

  ImageWriterConfig cfg = config;
  if (cfg.image_shape.x <= 0 || cfg.image_shape.y <= 0)
  {
    cfg.image_shape = array.shape;
  }
  if (cfg.tile_shape.x <= 0 || cfg.tile_shape.y <= 0)
  {
    cfg.tile_shape = array.shape;
  }

  std::unique_ptr<ImageWriter> writer = ImageWriter::create(fname, cfg);
  if (!writer || !writer->open(fname, cfg))
  {
    hmap::log::error("export_image failed to open writer for {}", fname);
    return false;
  }

  bool ok = writer->write_chunk(0, 0, array);
  writer->close();
  return ok;
}

bool export_image(const VirtualArray      &va,
                  const std::string       &fname,
                  const ImageWriterConfig &config)
{
  ComputeMode cm;
  cm.mode = ForEachMode::VA_DISTRIBUTED;
  return export_image(va, fname, config, cm);
}

bool export_image(const VirtualArray      &va,
                  const std::string       &fname,
                  const ImageWriterConfig &config,
                  const ComputeMode       &cm)
{
  ImageWriterConfig cfg = config;
  if (cfg.image_shape.x <= 0 || cfg.image_shape.y <= 0)
  {
    cfg.image_shape = va.shape;
  }
  if (cfg.tile_shape.x <= 0 || cfg.tile_shape.y <= 0)
  {
    cfg.tile_shape = va.tile_shape;
  }

  std::unique_ptr<ImageWriter> writer = ImageWriter::create(fname, cfg);
  if (!writer || !writer->open(fname, cfg))
  {
    hmap::log::error("export_image failed to open writer for {}", fname);
    return false;
  }

  return export_virtual_array(va, *writer, cm);
}

bool export_virtual_array(const VirtualArray &va, ImageWriter &writer)
{
  ComputeMode cm;
  cm.mode = ForEachMode::VA_DISTRIBUTED;
  return export_virtual_array(va, writer, cm);
}

bool export_virtual_array(const VirtualArray &va,
                          ImageWriter        &writer,
                          const ComputeMode  &cm)
{
  if (!writer.is_open())
  {
    hmap::log::error("export_virtual_array called with closed writer");
    return false;
  }

  for_each_tile(
      va,
      [&](const Array &tile, const TileRegion &region)
      { writer.write_tile(region, tile); },
      cm);

  writer.close();
  return true;
}

} // namespace hmap
