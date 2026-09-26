/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>

#include "highmap/export/opencv_writer.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

// =====================================
// Implementation details
// =====================================

struct OpenCVWriter::Impl
{
  std::string       filepath;
  ImageWriterConfig config;
  cv::Mat           mat;
  std::mutex        mutex;
  bool              is_open = false;
};

// =====================================
// OpenCVWriter class
// =====================================

OpenCVWriter::OpenCVWriter() : p_impl(std::make_unique<Impl>())
{
}

OpenCVWriter::OpenCVWriter(const std::string       &filepath,
                           const ImageWriterConfig &config)
    : p_impl(std::make_unique<Impl>())
{
  this->open(filepath, config);
}

OpenCVWriter::~OpenCVWriter()
{
  this->close();
}

void OpenCVWriter::close()
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (this->p_impl->is_open && !this->p_impl->mat.empty())
  {
    std::filesystem::path p(this->p_impl->filepath);
    if (p.has_parent_path())
    {
      std::filesystem::create_directories(p.parent_path());
    }

    cv::imwrite(this->p_impl->filepath, this->p_impl->mat);
    this->p_impl->mat.release();
    this->p_impl->is_open = false;
  }
}

const ImageWriterConfig &OpenCVWriter::config() const
{
  return this->p_impl->config;
}

const std::string &OpenCVWriter::filepath() const
{
  return this->p_impl->filepath;
}

bool OpenCVWriter::is_open() const
{
  return this->p_impl->is_open;
}

bool OpenCVWriter::open(const std::string       &filepath,
                        const ImageWriterConfig &config)
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (this->p_impl->is_open && !this->p_impl->mat.empty())
  {
    this->p_impl->mat.release();
    this->p_impl->is_open = false;
  }

  this->p_impl->filepath = filepath;
  this->p_impl->config = config;

  if (config.image_shape.x <= 0 || config.image_shape.y <= 0) return false;

  int cv_depth = CV_32F;
  if (config.data_type == ImageDataType::UINT8)
  {
    cv_depth = CV_8U;
  }
  else if (config.data_type == ImageDataType::UINT16)
  {
    cv_depth = CV_16U;
  }

  int cv_type = CV_MAKETYPE(cv_depth, config.channels);
  this->p_impl->mat = cv::Mat::zeros(config.image_shape.y,
                                     config.image_shape.x,
                                     cv_type);
  this->p_impl->is_open = true;
  return true;
}

bool OpenCVWriter::write_chunk(int offset_x, int offset_y, const Array &data)
{
  std::vector<Array> ch = {data};
  return this->write_chunk(offset_x, offset_y, ch);
}

bool OpenCVWriter::write_chunk(int                       offset_x,
                               int                       offset_y,
                               const std::vector<Array> &channels)
{
  std::lock_guard<std::mutex> lock(this->p_impl->mutex);

  if (!this->p_impl->is_open || channels.empty()) return false;

  const int src_w = channels[0].shape.x;
  const int src_h = channels[0].shape.y;
  const int nch = this->p_impl->config.channels;

  const bool  remap_val = this->p_impl->config.enable_remapping;
  const float vmin = this->p_impl->config.remap_min;
  const float vmax = this->p_impl->config.remap_max;
  const float vrange = (vmax > vmin) ? (vmax - vmin) : 1.f;

  for (int j = 0; j < src_h; ++j)
  {
    int dst_y = offset_y + j;
    if (dst_y < 0 || dst_y >= this->p_impl->mat.rows) continue;

    for (int i = 0; i < src_w; ++i)
    {
      int dst_x = offset_x + i;
      if (dst_x < 0 || dst_x >= this->p_impl->mat.cols) continue;

      if (this->p_impl->config.data_type == ImageDataType::FLOAT32)
      {
        if (nch == 1)
        {
          float val = channels[0](i, j);
          if (remap_val) val = (val - vmin) / vrange;
          this->p_impl->mat.at<float>(dst_y, dst_x) = val;
        }
        else
        {
          cv::Vec3f &pixel = this->p_impl->mat.at<cv::Vec3f>(dst_y, dst_x);
          for (int c = 0; c < 3 && c < static_cast<int>(channels.size()); ++c)
          {
            float val = channels[c](i, j);
            if (remap_val) val = (val - vmin) / vrange;
            pixel[c] = val;
          }
        }
      }
      else if (this->p_impl->config.data_type == ImageDataType::UINT16)
      {
        if (nch == 1)
        {
          float val = channels[0](i, j);
          if (remap_val) val = (val - vmin) / vrange;
          val = std::clamp(val, 0.f, 1.f) * 65535.f;
          this->p_impl->mat.at<uint16_t>(dst_y, dst_x) = static_cast<uint16_t>(
              std::round(val));
        }
        else
        {
          cv::Vec3w &pixel = this->p_impl->mat.at<cv::Vec3w>(dst_y, dst_x);
          for (int c = 0; c < 3 && c < static_cast<int>(channels.size()); ++c)
          {
            float val = channels[c](i, j);
            if (remap_val) val = (val - vmin) / vrange;
            val = std::clamp(val, 0.f, 1.f) * 65535.f;
            pixel[c] = static_cast<uint16_t>(std::round(val));
          }
        }
      }
      else if (this->p_impl->config.data_type == ImageDataType::UINT8)
      {
        if (nch == 1)
        {
          float val = channels[0](i, j);
          if (remap_val) val = (val - vmin) / vrange;
          val = std::clamp(val, 0.f, 1.f) * 255.f;
          this->p_impl->mat.at<uint8_t>(dst_y, dst_x) = static_cast<uint8_t>(
              std::round(val));
        }
        else
        {
          cv::Vec3b &pixel = this->p_impl->mat.at<cv::Vec3b>(dst_y, dst_x);
          for (int c = 0; c < 3 && c < static_cast<int>(channels.size()); ++c)
          {
            float val = channels[c](i, j);
            if (remap_val) val = (val - vmin) / vrange;
            val = std::clamp(val, 0.f, 1.f) * 255.f;
            pixel[c] = static_cast<uint8_t>(std::round(val));
          }
        }
      }
    }
  }

  return true;
}

bool OpenCVWriter::write_tile(const TileRegion &region, const Array &data)
{
  return ImageWriter::write_tile(region, data);
}

bool OpenCVWriter::write_tile(const TileRegion         &region,
                              const std::vector<Array> &channels)
{
  return ImageWriter::write_tile(region, channels);
}

} // namespace hmap
