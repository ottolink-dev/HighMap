/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "highmap/export/bigtiff_writer.hpp"
#include "highmap/logger.hpp"

#include <tiffio.h>

namespace hmap
{

// =====================================
// Implementation details
// =====================================

struct BigTiffWriter::Impl
{
  std::string       filepath;
  ImageWriterConfig config;
  TIFF             *tif = nullptr;
  std::mutex        mutex;
  bool              is_open = false;
};

// =====================================
// BigTiffWriter class
// =====================================

BigTiffWriter::BigTiffWriter() : p_impl(std::make_unique<Impl>())
{
}

BigTiffWriter::BigTiffWriter(const std::string       &filepath,
                             const ImageWriterConfig &config)
    : p_impl(std::make_unique<Impl>())
{
  this->open(filepath, config);
}

BigTiffWriter::~BigTiffWriter()
{
  this->close();
}

void BigTiffWriter::close()
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (this->p_impl->is_open && this->p_impl->tif != nullptr)
  {
    TIFFClose(this->p_impl->tif);
    this->p_impl->tif = nullptr;
    this->p_impl->is_open = false;
  }
}

const ImageWriterConfig &BigTiffWriter::config() const
{
  return this->p_impl->config;
}

const std::string &BigTiffWriter::filepath() const
{
  return this->p_impl->filepath;
}

bool BigTiffWriter::is_open() const
{
  return this->p_impl->is_open;
}

bool BigTiffWriter::open(const std::string       &filepath,
                         const ImageWriterConfig &config)
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (this->p_impl->is_open && this->p_impl->tif != nullptr)
  {
    TIFFClose(this->p_impl->tif);
    this->p_impl->tif = nullptr;
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

  // open in BigTIFF mode ("w8")
  TIFF *tif = TIFFOpen(filepath.c_str(), "w8");
  if (!tif)
  {
    hmap::log::error("BigTiffWriter failed to open file: {}", filepath);
    return false;
  }

  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, config.image_shape.x);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, config.image_shape.y);
  TIFFSetField(tif, TIFFTAG_TILEWIDTH, config.tile_shape.x);
  TIFFSetField(tif, TIFFTAG_TILELENGTH, config.tile_shape.y);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, config.channels);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

  if (config.channels >= 3)
  {
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  }
  else
  {
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  }

  // bit depth and format
  if (config.data_type == ImageDataType::FLOAT32)
  {
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  }
  else if (config.data_type == ImageDataType::UINT16)
  {
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
  }
  else if (config.data_type == ImageDataType::UINT8)
  {
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
  }
  else
  {
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  }

  // compression
  std::string comp = config.compression;
  std::transform(comp.begin(), comp.end(), comp.begin(), ::tolower);

  if (comp == "lzw")
  {
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_LZW);
  }
  else if (comp == "deflate" || comp == "zip" || comp == "zlib")
  {
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
  }
  else if (comp == "packbits")
  {
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
  }
  else
  {
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
  }

  this->p_impl->tif = tif;
  this->p_impl->is_open = true;
  return true;
}

bool BigTiffWriter::write_chunk(int offset_x, int offset_y, const Array &data)
{
  std::vector<Array> ch = {data};
  return this->write_chunk(offset_x, offset_y, ch);
}

bool BigTiffWriter::write_chunk(int                       offset_x,
                                int                       offset_y,
                                const std::vector<Array> &channels)
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (!this->p_impl->is_open || !this->p_impl->tif) return false;
  if (channels.empty()) return false;

  const int tw = this->p_impl->config.tile_shape.x;
  const int th = this->p_impl->config.tile_shape.y;
  const int nch = this->p_impl->config.channels;

  const int src_w = channels[0].shape.x;
  const int src_h = channels[0].shape.y;

  const bool  remap_val = this->p_impl->config.enable_remapping;
  const float vmin = this->p_impl->config.remap_min;
  const float vmax = this->p_impl->config.remap_max;
  const float vrange = (vmax > vmin) ? (vmax - vmin) : 1.f;

  if (this->p_impl->config.data_type == ImageDataType::FLOAT32)
  {
    std::vector<float> tile_buf(tw * th * nch, 0.f);

    for (int j = 0; j < src_h && j < th; ++j)
    {
      for (int i = 0; i < src_w && i < tw; ++i)
      {
        size_t dst_idx = (j * tw + i) * nch;
        for (int c = 0; c < nch && c < static_cast<int>(channels.size()); ++c)
        {
          float val = channels[c](i, j);
          if (remap_val) val = (val - vmin) / vrange;
          tile_buf[dst_idx + c] = val;
        }
      }
    }

    tmsize_t written = TIFFWriteTile(this->p_impl->tif,
                                     tile_buf.data(),
                                     offset_x,
                                     offset_y,
                                     0,
                                     0);
    return (written >= 0);
  }
  else if (this->p_impl->config.data_type == ImageDataType::UINT16)
  {
    std::vector<uint16_t> tile_buf(tw * th * nch, 0);

    for (int j = 0; j < src_h && j < th; ++j)
    {
      for (int i = 0; i < src_w && i < tw; ++i)
      {
        size_t dst_idx = (j * tw + i) * nch;
        for (int c = 0; c < nch && c < static_cast<int>(channels.size()); ++c)
        {
          float val = channels[c](i, j);
          if (remap_val)
          {
            val = (val - vmin) / vrange;
          }
          val = std::clamp(val, 0.f, 1.f) * 65535.f;
          tile_buf[dst_idx + c] = static_cast<uint16_t>(std::round(val));
        }
      }
    }

    tmsize_t written = TIFFWriteTile(this->p_impl->tif,
                                     tile_buf.data(),
                                     offset_x,
                                     offset_y,
                                     0,
                                     0);
    return (written >= 0);
  }
  else if (this->p_impl->config.data_type == ImageDataType::UINT8)
  {
    std::vector<uint8_t> tile_buf(tw * th * nch, 0);

    for (int j = 0; j < src_h && j < th; ++j)
    {
      for (int i = 0; i < src_w && i < tw; ++i)
      {
        size_t dst_idx = (j * tw + i) * nch;
        for (int c = 0; c < nch && c < static_cast<int>(channels.size()); ++c)
        {
          float val = channels[c](i, j);
          if (remap_val)
          {
            val = (val - vmin) / vrange;
          }
          val = std::clamp(val, 0.f, 1.f) * 255.f;
          tile_buf[dst_idx + c] = static_cast<uint8_t>(std::round(val));
        }
      }
    }

    tmsize_t written = TIFFWriteTile(this->p_impl->tif,
                                     tile_buf.data(),
                                     offset_x,
                                     offset_y,
                                     0,
                                     0);
    return (written >= 0);
  }

  return false;
}

bool BigTiffWriter::write_tile(const TileRegion &region, const Array &data)
{
  return ImageWriter::write_tile(region, data);
}

} // namespace hmap
