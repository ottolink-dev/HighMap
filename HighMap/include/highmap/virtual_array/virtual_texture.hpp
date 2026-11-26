/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file virtual_array.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <memory>
#include <string>
#include <vector>

#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

#include "highmap/colorize.hpp"
#include "highmap/texture.hpp"
#include "highmap/virtual_array/virtual_array.hpp"
#include "highmap/virtual_array/virtual_texture_storage.hpp"

namespace hmap
{

class VirtualTexture
{
public:
  // --- Ctor

  VirtualTexture() = default;

  VirtualTexture(glm::ivec2                   shape,
                 glm::ivec2                   tile_shape,
                 int                          halo,
                 int                          channels,
                 std::unique_ptr<TileStorage> storage_proto);

  VirtualTexture(glm::ivec2  shape,
                 glm::ivec2  tile_shape,
                 int         halo,
                 int         channels,
                 StorageMode storage_mode);

  // --- Copy

  void copy_from(VirtualTexture &src, const ComputeMode &cm);

  // --- Channels data

  int                               channels() const;
  VirtualArray                     &channel(int c);
  const VirtualArray               &channel(int c) const;
  std::vector<VirtualArray *>       channels_ptr();
  std::vector<const VirtualArray *> channels_ptr() const;
  std::vector<VirtualArray>        &get_arrays();

  void fill(float value, const ComputeMode &cm);
  void fill(int c, float value, const ComputeMode &cm);

  void sync_overlap_buffers(SyncOperation op = SyncOperation::SmoothBlend);

  // --- Converters

  void from_arrays(const std::vector<const Array *> &p_arrays,
                   const ComputeMode                &cm);

  std::vector<uint8_t> to_img_8bit(const glm::ivec2  &img_shape,
                                   const ComputeMode &c,
                                   bool               flip_y = false) const;

  void to_png(const glm::ivec2  &array_shape,
              const std::string &fname,
              const ComputeMode &cm,
              int                depth = CV_8U) const;

  void to_png(const std::string &fname,
              const ComputeMode &cm,
              int                depth = CV_8U) const;

  std::vector<float> to_raw(const ComputeMode &cm, bool flip_y = false);

  Texture to_texture(const glm::ivec2 &img_shape, const ComputeMode &cm) const;

  void   trim_storage() const;
  size_t live_memory_bytes() const;

  // --- Members

  glm::ivec2 shape;
  glm::ivec2 tile_shape;
  int        halo;

private:
  std::vector<VirtualArray> arrays;
};

// functions

VirtualTexture convert_texture_channels(const VirtualTexture &src,
                                        int                   dst_channels,
                                        float                 fill_value,
                                        const ComputeMode    &cm);

template <typename Func>
void for_each_tile(VirtualTexture &tex, Func &&func, const ComputeMode &cm);

template <typename Func>
void for_each_tile(const VirtualTexture &tex,
                   Func                &&func,
                   const ComputeMode    &cm);

template <typename Func>
void for_each_pixel(VirtualTexture &tex, Func &&func, const ComputeMode &cm);

#include "highmap/virtual_array/virtual_texture.inl"

} // namespace hmap
