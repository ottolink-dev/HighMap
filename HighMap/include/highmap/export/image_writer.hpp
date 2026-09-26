/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file image_writer.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for generic image writer interface.
 *
 * @copyright Copyright (c) 2026 Otto Link
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/virtual_array/tile_region.hpp"

namespace hmap
{

// =====================================
// Data types and configuration
// =====================================

enum class ImageDataType : int
{
  FLOAT32,
  FLOAT16,
  UINT16,
  UINT8
};

struct ImageWriterConfig
{
  glm::ivec2    image_shape = {0, 0};
  glm::ivec2    tile_shape = {512, 512};
  ImageDataType data_type = ImageDataType::FLOAT32;
  int           channels = 1;
  std::string   compression = "default";
  float         remap_min = 0.f;
  float         remap_max = 1.f;
  bool          enable_remapping = false;
};

// =====================================
// Abstract class
// =====================================

class ImageWriter
{
public:
  virtual ~ImageWriter() = default;

  // --- Lifecycle

  virtual bool open(const std::string       &filepath,
                    const ImageWriterConfig &config) = 0;
  virtual void close() = 0;
  virtual bool is_open() const = 0;

  // --- Chunk / tile writing

  virtual bool write_chunk(int offset_x, int offset_y, const Array &data) = 0;
  virtual bool write_chunk(int                       offset_x,
                           int                       offset_y,
                           const std::vector<Array> &channels);
  virtual bool write_tile(const TileRegion &region, const Array &data);

  // --- Properties

  virtual const ImageWriterConfig &config() const = 0;
  virtual const std::string       &filepath() const = 0;

  // --- Factory

  static std::unique_ptr<ImageWriter> create(const std::string       &filepath,
                                             const ImageWriterConfig &config);
};

} // namespace hmap
