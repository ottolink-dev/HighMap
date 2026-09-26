/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "highmap/export/openexr_writer.hpp"
#include "highmap/logger.hpp"

#include <Imath/ImathBox.h>
#include <Imath/half.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfTileDescription.h>
#include <OpenEXR/ImfTiledOutputFile.h>

namespace hmap
{

// =====================================
// Implementation details
// =====================================

struct OpenEXRWriter::Impl
{
  std::string                           filepath;
  ImageWriterConfig                     config;
  std::unique_ptr<Imf::TiledOutputFile> file;
  std::mutex                            mutex;
  bool                                  is_open = false;
};

// =====================================
// OpenEXRWriter class
// =====================================

OpenEXRWriter::OpenEXRWriter() : p_impl(std::make_unique<Impl>())
{
}

OpenEXRWriter::OpenEXRWriter(const std::string       &filepath,
                             const ImageWriterConfig &config)
    : p_impl(std::make_unique<Impl>())
{
  this->open(filepath, config);
}

OpenEXRWriter::~OpenEXRWriter()
{
  this->close();
}

void OpenEXRWriter::close()
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (this->p_impl->is_open)
  {
    this->p_impl->file.reset();
    this->p_impl->is_open = false;
  }
}

const ImageWriterConfig &OpenEXRWriter::config() const
{
  return this->p_impl->config;
}

const std::string &OpenEXRWriter::filepath() const
{
  return this->p_impl->filepath;
}

bool OpenEXRWriter::is_open() const
{
  return this->p_impl->is_open;
}

bool OpenEXRWriter::open(const std::string       &filepath,
                         const ImageWriterConfig &config)
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (this->p_impl->is_open)
  {
    this->p_impl->file.reset();
    this->p_impl->is_open = false;
  }

  this->p_impl->filepath = filepath;
  this->p_impl->config = config;

  if (config.image_shape.x <= 0 || config.image_shape.y <= 0) return false;
  if (config.tile_shape.x <= 0 || config.tile_shape.y <= 0) return false;

  // create parent directory if needed
  std::filesystem::path p(filepath);
  if (p.has_parent_path())
  {
    std::filesystem::create_directories(p.parent_path());
  }

  try
  {
    Imath::Box2i display_window(
        Imath::V2i(0, 0),
        Imath::V2i(config.image_shape.x - 1, config.image_shape.y - 1));
    Imath::Box2i data_window(
        Imath::V2i(0, 0),
        Imath::V2i(config.image_shape.x - 1, config.image_shape.y - 1));

    Imf::Header header(display_window, data_window);
    header.lineOrder() = Imf::RANDOM_Y;

    header.setTileDescription(Imf::TileDescription(config.tile_shape.x,
                                                   config.tile_shape.y,
                                                   Imf::ONE_LEVEL));

    Imf::PixelType pixel_type = (config.data_type == ImageDataType::FLOAT16)
                                    ? Imf::HALF
                                    : Imf::FLOAT;

    if (config.channels == 1)
    {
      header.channels().insert("Y", Imf::Channel(pixel_type));
    }
    else if (config.channels == 3)
    {
      header.channels().insert("R", Imf::Channel(pixel_type));
      header.channels().insert("G", Imf::Channel(pixel_type));
      header.channels().insert("B", Imf::Channel(pixel_type));
    }
    else if (config.channels == 4)
    {
      header.channels().insert("R", Imf::Channel(pixel_type));
      header.channels().insert("G", Imf::Channel(pixel_type));
      header.channels().insert("B", Imf::Channel(pixel_type));
      header.channels().insert("A", Imf::Channel(pixel_type));
    }
    else
    {
      return false;
    }

    this->p_impl->file = std::make_unique<Imf::TiledOutputFile>(
        filepath.c_str(),
        header);
    this->p_impl->is_open = true;
    return true;
  }
  catch (const std::exception &e)
  {
    hmap::log::error("OpenEXRWriter failed to open {}: {}", filepath, e.what());
    this->p_impl->is_open = false;
    return false;
  }
}

bool OpenEXRWriter::write_chunk(int offset_x, int offset_y, const Array &data)
{
  std::vector<Array> ch = {data};
  return this->write_chunk(offset_x, offset_y, ch);
}

bool OpenEXRWriter::write_chunk(int                       offset_x,
                                int                       offset_y,
                                const std::vector<Array> &channels)
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (!this->p_impl->is_open || !this->p_impl->file) return false;
  if (channels.empty()) return false;

  const int tw = this->p_impl->config.tile_shape.x;
  const int th = this->p_impl->config.tile_shape.y;

  if (tw <= 0 || th <= 0) return false;

  int tile_x = offset_x / tw;
  int tile_y = offset_y / th;

  try
  {
    Imf::FrameBuffer               fb;
    std::vector<std::vector<half>> half_buffers;

    const size_t                          num_channels = channels.size();
    static const std::vector<std::string> names_1 = {"Y"};
    static const std::vector<std::string> names_3 = {"R", "G", "B"};
    static const std::vector<std::string> names_4 = {"R", "G", "B", "A"};

    const std::vector<std::string> *p_names = &names_1;
    if (this->p_impl->config.channels == 3) p_names = &names_3;
    if (this->p_impl->config.channels == 4) p_names = &names_4;

    const bool is_half = (this->p_impl->config.data_type ==
                          ImageDataType::FLOAT16);

    const bool  remap_val = this->p_impl->config.enable_remapping;
    const float vmin = this->p_impl->config.remap_min;
    const float vmax = this->p_impl->config.remap_max;
    const float vrange = (vmax > vmin) ? (vmax - vmin) : 1.f;

    std::vector<std::vector<float>> float_remapped_buffers;

    for (size_t c = 0; c < p_names->size(); ++c)
    {
      const Array       &arr = (c < num_channels) ? channels[c] : channels[0];
      const std::string &chan_name = (*p_names)[c];

      if (is_half)
      {
        half_buffers.emplace_back(arr.vector.size());
        auto &hbuf = half_buffers.back();
        for (size_t k = 0; k < arr.vector.size(); ++k)
        {
          float v = arr.vector[k];
          if (remap_val) v = (v - vmin) / vrange;
          hbuf[k] = half(v);
        }

        char *base = reinterpret_cast<char *>(hbuf.data()) -
                     (static_cast<ptrdiff_t>(offset_x) * sizeof(half) +
                      static_cast<ptrdiff_t>(offset_y) * arr.shape.x *
                          sizeof(half));

        fb.insert(chan_name,
                  Imf::Slice(Imf::HALF,
                             base,
                             sizeof(half),
                             arr.shape.x * sizeof(half)));
      }
      else
      {
        const float *data_ptr = arr.vector.data();
        if (remap_val)
        {
          float_remapped_buffers.emplace_back(arr.vector.size());
          auto &fbuf = float_remapped_buffers.back();
          for (size_t k = 0; k < arr.vector.size(); ++k)
          {
            fbuf[k] = (arr.vector[k] - vmin) / vrange;
          }
          data_ptr = fbuf.data();
        }

        char *base = reinterpret_cast<char *>(const_cast<float *>(data_ptr)) -
                     (static_cast<ptrdiff_t>(offset_x) * sizeof(float) +
                      static_cast<ptrdiff_t>(offset_y) * arr.shape.x *
                          sizeof(float));

        fb.insert(chan_name,
                  Imf::Slice(Imf::FLOAT,
                             base,
                             sizeof(float),
                             arr.shape.x * sizeof(float)));
      }
    }

    this->p_impl->file->setFrameBuffer(fb);
    this->p_impl->file->writeTile(tile_x, tile_y);
    return true;
  }
  catch (const std::exception &e)
  {
    hmap::log::error("OpenEXRWriter failed to write tile ({}, {}): {}",
                     tile_x,
                     tile_y,
                     e.what());
    return false;
  }
}

bool OpenEXRWriter::write_tile(const TileRegion &region, const Array &data)
{
  return ImageWriter::write_tile(region, data);
}

bool OpenEXRWriter::write_tile(const TileRegion         &region,
                               const std::vector<Array> &channels)
{
  return ImageWriter::write_tile(region, channels);
}

} // namespace hmap
