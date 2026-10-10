/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/texture.hpp"

#include <format>

namespace hmap
{

Texture::Texture(glm::ivec2 shape, int num_channels)
    : shape(shape), channels(num_channels, Array(shape))
{
}

Texture::Texture(glm::ivec2 shape, int num_channels, float fill_value)
    : shape(shape), channels(num_channels, Array(shape, fill_value))
{
}

Texture::Texture(const Array &array) : shape(array.shape), channels({array})
{
}

Texture::Texture(const Array &r, const Array &g, const Array &b)
    : shape(r.shape), channels({r, g, b})
{
  if (g.shape != shape || b.shape != shape)
  {
    throw std::runtime_error("Texture: Input arrays must have matching shape");
  }
}

Texture::Texture(const Array &r, const Array &g, const Array &b, const Array &a)
    : shape(r.shape), channels({r, g, b, a})
{
  if (g.shape != shape || b.shape != shape || a.shape != shape)
  {
    throw std::runtime_error("Texture: Input arrays must have matching shape");
  }
}

Texture::Texture(const std::vector<Array> &channels) : channels(channels)
{
  if (!channels.empty())
  {
    shape = channels[0].shape;
    for (const auto &ch : channels)
    {
      if (ch.shape != shape)
      {
        throw std::runtime_error(
            "Texture: Input arrays must have matching shape");
      }
    }
  }
}

Texture::Texture(const std::string &fname, bool flip_j)
{
  cv::Mat mat = cv::imread(fname, cv::IMREAD_COLOR);

  if (mat.data == nullptr)
  {
    hmap::log::error("error while reading the image file: {}", fname);
    this->shape = glm::ivec2(0, 0);
    this->channels.clear();
    return;
  }

  cv::cvtColor(mat, mat, cv::COLOR_BGR2RGB);
  mat.convertTo(mat, CV_32FC3, 1.f / 255.f);

  this->shape = glm::ivec2(mat.cols, mat.rows);
  this->channels.clear();
  // 4 channels (RGBA)
  this->channels.resize(4, Array(this->shape));

  // fill texture
  for (int j = 0; j < this->shape.y; j++)
  {
    int jj = flip_j ? (this->shape.y - 1 - j) : j;
    for (int i = 0; i < this->shape.x; i++)
    {
      cv::Vec3f pixel = mat.at<cv::Vec3f>(j, i);

      // assign RGB values to planar arrays
      this->channels[0](i, jj) = pixel[0]; // red
      this->channels[1](i, jj) = pixel[1]; // green
      this->channels[2](i, jj) = pixel[2]; // blue
      this->channels[3](i, jj) = 1.f;      // alpha
    }
  }
}

glm::vec3 Texture::get_pixel3(int i, int j) const
{
  if (num_channels() < 3)
  {
    throw std::runtime_error(
        "Texture::get_pixel3: Texture has fewer than 3 channels");
  }
  return glm::vec3(channels[0](i, j), channels[1](i, j), channels[2](i, j));
}

glm::vec4 Texture::get_pixel4(int i, int j) const
{
  if (num_channels() < 4)
  {
    throw std::runtime_error(
        "Texture::get_pixel4: Texture has fewer than 4 channels");
  }
  return glm::vec4(channels[0](i, j),
                   channels[1](i, j),
                   channels[2](i, j),
                   channels[3](i, j));
}

void Texture::set_pixel(int i, int j, const glm::vec3 &color)
{
  if (num_channels() < 3)
  {
    throw std::runtime_error(
        "Texture::set_pixel: Texture has fewer than 3 channels");
  }
  channels[0](i, j) = color.x;
  channels[1](i, j) = color.y;
  channels[2](i, j) = color.z;
}

void Texture::set_pixel(int i, int j, const glm::vec4 &color)
{
  if (num_channels() < 4)
  {
    throw std::runtime_error(
        "Texture::set_pixel: Texture has fewer than 4 channels");
  }
  channels[0](i, j) = color.x;
  channels[1](i, j) = color.y;
  channels[2](i, j) = color.z;
  channels[3](i, j) = color.w;
}

float Texture::max() const
{
  if (channels.empty()) return 0.f;
  float max_val = channels[0].max();
  for (size_t c = 1; c < channels.size(); ++c)
  {
    max_val = std::max(max_val, channels[c].max());
  }
  return max_val;
}

float Texture::min() const
{
  if (channels.empty()) return 0.f;
  float min_val = channels[0].min();
  for (size_t c = 1; c < channels.size(); ++c)
  {
    min_val = std::min(min_val, channels[c].min());
  }
  return min_val;
}

void Texture::remap(float vmin, float vmax)
{
  if (!validate_non_empty(*this)) return;

  float min_val = this->min();
  float max_val = this->max();

  if (min_val != max_val)
  {
    for (auto &ch : channels)
    {
      for (auto &v : ch.vector)
      {
        v = (v - min_val) / (max_val - min_val) * (vmax - vmin) + vmin;
      }
    }
  }
  else
  {
    for (auto &ch : channels)
    {
      std::fill(ch.vector.begin(), ch.vector.end(), vmin);
    }
  }
}

Texture Texture::resample_to_shape(glm::ivec2 new_shape_xy) const
{
  if (!validate_non_empty(*this) || !validate_shape(new_shape_xy))
    return Texture();

  Texture out;
  out.shape = new_shape_xy;
  out.channels.reserve(this->channels.size());
  for (const auto &ch : this->channels)
  {
    out.channels.push_back(ch.resample_to_shape(new_shape_xy));
  }
  return out;
}

cv::Mat Texture::to_cv_mat() const
{
  if (!validate_non_empty(*this)) return cv::Mat();

  int nch = num_channels();

  std::vector<cv::Mat> cv_channels;
  cv_channels.reserve(nch);
  for (int c = 0; c < nch; ++c)
  {
    // Need to cast away const for cv::Mat constructor, but cv::merge copies the
    // data.
    cv_channels.push_back(
        cv::Mat(shape.y,
                shape.x,
                CV_32FC1,
                const_cast<float *>(channels[c].vector.data())));
  }

  cv::Mat merged;
  cv::merge(cv_channels, merged);

  if (nch == 3)
  {
    cv::cvtColor(merged, merged, cv::COLOR_RGB2BGR);
  }
  else if (nch == 4)
  {
    cv::cvtColor(merged, merged, cv::COLOR_RGBA2BGRA);
  }

  return merged;
}

void Texture::to_png(const std::string &fname, int depth) const
{
  if (!validate_non_empty(*this)) return;

  cv::Mat mat = to_cv_mat();
  int     scale_factor = (depth == CV_8U) ? 255 : 65535;
  mat.convertTo(mat, depth, scale_factor);
  cv::flip(mat, mat, 0); // up-down
  cv::imwrite(fname, mat);
}

std::vector<uint8_t> Texture::to_img_8bit(bool flip_y) const
{
  if (!validate_non_empty(*this)) return {};

  std::vector<uint8_t> vec;
  int                  nch = num_channels();
  vec.reserve(shape.x * shape.y * nch);

  if (flip_y)
  {
    for (int j = shape.y - 1; j >= 0; j--)
      for (int i = 0; i < shape.x; i++)
        for (int c = 0; c < nch; c++)
          vec.push_back(static_cast<uint8_t>(
              std::clamp(255.f * channels[c](i, j), 0.f, 255.f)));
  }
  else
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; i++)
        for (int c = 0; c < nch; c++)
          vec.push_back(static_cast<uint8_t>(
              std::clamp(255.f * channels[c](i, j), 0.f, 255.f)));
  }

  return vec;
}

} // namespace hmap
