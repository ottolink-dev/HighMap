/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file opencv_writer.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief OpenCV-based image writer implementing ImageWriter interface.
 *
 * @copyright Copyright (c) 2026 Otto Link
 */
#pragma once

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "highmap/export/image_writer.hpp"

namespace hmap
{

class OpenCVWriter : public ImageWriter
{
public:
  OpenCVWriter();
  explicit OpenCVWriter(const std::string       &filepath,
                        const ImageWriterConfig &config);
  ~OpenCVWriter() override;

  // --- Lifecycle

  bool open(const std::string       &filepath,
            const ImageWriterConfig &config) override;
  void close() override;
  bool is_open() const override;

  // --- Chunk / tile writing

  bool write_chunk(int offset_x, int offset_y, const Array &data) override;
  bool write_chunk(int                       offset_x,
                   int                       offset_y,
                   const std::vector<Array> &channels) override;
  bool write_tile(const TileRegion &region, const Array &data) override;

  // --- Properties

  const ImageWriterConfig &config() const override;
  const std::string       &filepath() const override;

private:
  struct Impl;
  std::unique_ptr<Impl> p_impl;
};

} // namespace hmap
