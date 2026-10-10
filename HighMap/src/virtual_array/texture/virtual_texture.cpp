/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/texture.hpp"
#include "highmap/virtual_array/tile_storage.hpp"
#include "highmap/virtual_array/virtual_array.hpp"
#include "highmap/virtual_array/virtual_texture.hpp"

#include <format>

namespace hmap
{

VirtualTexture::VirtualTexture(glm::ivec2                   shape,
                               glm::vec4                    bbox,
                               glm::ivec2                   tile_shape,
                               int                          halo,
                               int                          channels,
                               std::unique_ptr<TileStorage> storage_proto)
    : shape(shape), bbox(bbox), tile_shape(tile_shape), halo(halo)
{
  arrays.reserve(channels);

  for (int c = 0; c < channels; ++c)
  {
    arrays.emplace_back(shape, bbox, tile_shape, halo, storage_proto->clone());
  }
}

VirtualTexture::VirtualTexture(glm::ivec2  shape,
                               glm::vec4   bbox,
                               glm::ivec2  tile_shape,
                               int         halo,
                               int         channels,
                               StorageMode storage_mode)
    : shape(shape), bbox(bbox), tile_shape(tile_shape), halo(halo)
{
  arrays.reserve(channels);

  for (int c = 0; c < channels; ++c)
  {
    arrays.emplace_back(shape, bbox, tile_shape, halo, storage_mode);
  }
}

VirtualTexture::VirtualTexture(glm::ivec2  shape,
                               glm::ivec2  tile_shape,
                               int         halo,
                               int         channels,
                               StorageMode storage_mode)
    : shape(shape),
      bbox({0.f, 1.f, 0.f, 1.f}),
      tile_shape(tile_shape),
      halo(halo)
{
  arrays.reserve(channels);

  for (int c = 0; c < channels; ++c)
  {
    arrays.emplace_back(shape, bbox, tile_shape, halo, storage_mode);
  }
}

int VirtualTexture::channels() const
{
  return int(arrays.size());
}

VirtualArray &VirtualTexture::channel(int c)
{
  return arrays[c];
}
const VirtualArray &VirtualTexture::channel(int c) const
{
  return arrays[c];
}

std::vector<VirtualArray *> VirtualTexture::channels_ptr()
{
  std::vector<VirtualArray *> out;
  out.reserve(arrays.size());

  for (auto &a : arrays)
    out.push_back(&a);

  return out;
}

std::vector<const VirtualArray *> VirtualTexture::channels_ptr() const
{
  std::vector<const VirtualArray *> out;
  out.reserve(arrays.size());

  for (const auto &a : arrays)
    out.push_back(&a);

  return out;
}

void VirtualTexture::copy_from(VirtualTexture &src, const ComputeMode &cm)
{
  if (this == &src) return;

  shape = src.shape;
  bbox = src.bbox;
  tile_shape = src.tile_shape;
  halo = src.halo;

  this->arrays.resize(src.channels());

  for (int c = 0; c < this->channels(); ++c)
    this->channel(c).copy_from(src.channel(c), cm);
}

void VirtualTexture::fill(float value, const ComputeMode &cm)
{
  for (auto &va : this->get_arrays())
    va.fill(value, cm);
}

void VirtualTexture::fill(int c, float value, const ComputeMode &cm)
{
  if (c > this->channels() - 1) return;

  this->channel(c).fill(value, cm);
}

void VirtualTexture::from_arrays(const std::vector<const Array *> &p_arrays,
                                 const ComputeMode                &cm)
{
  if (!validate_non_empty(p_arrays, "Array pointer list")) return;

  const int nch = this->channels();

  if (nch != int(p_arrays.size()))
  {
    hmap::log::error("size mismatch between arrays and channels nb ({} != {})",
                     nch,
                     p_arrays.size());
    return;
  }

  for (int c = 0; c < nch; ++c)
    if (p_arrays[c]) this->channel(c).from_array(*p_arrays[c], cm);
}

std::vector<VirtualArray> &VirtualTexture::get_arrays()
{
  return this->arrays;
}

std::vector<uint8_t> VirtualTexture::to_img_8bit(const glm::ivec2  &img_shape,
                                                 const ComputeMode &cm,
                                                 bool flip_y) const
{
  const int nch = this->channels();

  if (nch != 3 && nch != 4)
  {
    hmap::log::error("to_img_8bit supports only RGB or RGBA (got {})", nch);
    return {};
  }

  Texture t = this->to_texture(img_shape, cm);
  return t.to_img_8bit(flip_y);
}

void VirtualTexture::to_png(const glm::ivec2  &img_shape,
                            const std::string &fname,
                            const ComputeMode &cm,
                            int                depth) const
{
  const int nch = this->channels();

  if (nch != 3 && nch != 4)
  {
    hmap::log::error("PNG supports only RGB or RGBA (got {})", nch);
    return;
  }

  Texture t = this->to_texture(img_shape, cm);
  t.to_png(fname, depth);
}

void VirtualTexture::to_png(const std::string &fname,
                            const ComputeMode &cm,
                            int                depth) const
{
  this->to_png(this->shape, fname, cm, depth);
}

std::vector<float> VirtualTexture::to_raw(const ComputeMode &cm, bool flip_y)
{
  std::vector<float> raw;
  raw.reserve(this->shape.x * this->shape.y * this->channels());

  // brute force
  Texture t = this->to_texture(this->shape, cm);

  if (flip_y)
  {
    for (int j = 0; j < this->shape.y; ++j)
      for (int i = 0; i < this->shape.x; ++i)
        for (int c = 0; c < int(this->channels()); ++c)
          raw.push_back(t(i, this->shape.y - 1 - j, c));
  }
  else
  {
    for (int j = 0; j < this->shape.y; ++j)
      for (int i = 0; i < this->shape.x; ++i)
        for (int c = 0; c < int(this->channels()); ++c)
          raw.push_back(t(i, j, c));
  }

  return raw;
}

Texture VirtualTexture::to_texture(const glm::ivec2  &img_shape,
                                   const ComputeMode &cm) const
{
  if (!validate_shape(img_shape)) return Texture();

  const int nch = this->channels();

  if (nch < 1) throw std::runtime_error("VirtualTexture has no channels");

  Texture t(img_shape, nch);

  for (int c = 0; c < nch; ++c)
  {
    t[c] = this->channel(c).to_array(img_shape, cm);
  }

  return t;
}

void VirtualTexture::trim_storage() const
{
  for (const auto &va : this->arrays)
    va.trim_storage();
}

size_t VirtualTexture::live_memory_bytes() const
{
  size_t bytes = 0;
  for (const auto &va : this->arrays)
    bytes += va.live_memory_bytes();
  return bytes;
}

} // namespace hmap
