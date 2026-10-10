/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/tile_storage.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <format>

namespace hmap
{

VirtualArray::VirtualArray(glm::ivec2                   shape,
                           glm::vec4                    bbox,
                           glm::ivec2                   tile_shape,
                           int                          halo,
                           std::unique_ptr<TileStorage> storage)
    : shape(shape),
      bbox(bbox),
      tile_shape(tile_shape),
      halo(halo),
      storage(std::move(storage))
{
}

VirtualArray::VirtualArray(glm::ivec2  shape,
                           glm::vec4   bbox,
                           glm::ivec2  tile_shape,
                           int         halo,
                           StorageMode storage_mode)
    : shape(shape), bbox(bbox), tile_shape(tile_shape), halo(halo)
{
  this->storage = make_storage(this->shape, this->tile_shape, storage_mode);
}

VirtualArray::VirtualArray(glm::ivec2  shape,
                           glm::ivec2  tile_shape,
                           int         halo,
                           StorageMode storage_mode)
    : shape(shape),
      bbox({0.f, 1.f, 0.f, 1.f}),
      tile_shape(tile_shape),
      halo(halo)
{
  this->storage = make_storage(this->shape, this->tile_shape, storage_mode);
}

void VirtualArray::copy_from(const VirtualArray &src,
                             const ComputeMode  &cm,
                             bool                copy_src_data)
{
  if (this == &src) return;

  shape = src.shape;
  bbox = src.bbox;
  tile_shape = src.tile_shape;
  halo = src.halo;

  storage = src.storage->clone();

  if (copy_src_data) copy_data(src, *this, cm);
}

std::unique_ptr<VirtualArray> VirtualArray::clone(const ComputeMode &cm,
                                                  bool               deep_copy)
{
  // create empty VA with same geometry
  auto va = std::make_unique<VirtualArray>(
      this->shape,
      this->bbox,
      this->tile_shape,
      this->halo,
      StorageMode::VA_RAM // temporary, replaced below
  );

  // clone storage (empty)
  va->storage = this->storage->clone();

  // copy data if requested
  if (deep_copy) copy_data(*this, *va, cm);

  return va;
}

bool VirtualArray::empty() const
{
  return this->shape.x * this->shape.y == 0;
}

void VirtualArray::fill(float value, const ComputeMode &cm)
{
  for_each_tile(
      *this,
      [value](Array &tile, const TileRegion &) { tile = value; },
      cm);
}

void VirtualArray::from_array(const Array &array, const ComputeMode &cm)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [&array, this](Array &tile, const TileRegion &region)
  {
    const float rx = float(array.shape.x - 1) / float(this->shape.x - 1);
    const float ry = float(array.shape.y - 1) / float(this->shape.y - 1);

    glm::ivec2 ij0 = this->tile_region_global_indices(region);

    for (int j = 0; j < region.shape.y; ++j)
      for (int i = 0; i < region.shape.x; ++i)
      {
        // destination pixel center
        float xd = float(ij0.x + i) + 0.5f;
        float yd = float(ij0.y + j) + 0.5f;

        // source pixel center
        int ig = int(std::floor(xd * rx));
        int jg = int(std::floor(yd * ry));

        // clamp
        ig = std::clamp(ig, 0, array.shape.x - 1);
        jg = std::clamp(jg, 0, array.shape.y - 1);

        tile(i, j) = array(ig, jg);
      }
  };

  for_each_tile(*this, lambda, cm);
}

void VirtualArray::from_array_bilinear(const Array       &array,
                                       const ComputeMode &cm)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [&array, this](Array &tile, const TileRegion &region)
  {
    const float rx = float(array.shape.x - 1) / float(this->shape.x - 1);
    const float ry = float(array.shape.y - 1) / float(this->shape.y - 1);

    const glm::ivec2 ij0 = this->tile_region_global_indices(region);

    for (int j = 0; j < region.shape.y; ++j)
      for (int i = 0; i < region.shape.x; ++i)
      {
        // destination pixel center (global)
        const float xd = float(ij0.x + i) + 0.5f;
        const float yd = float(ij0.y + j) + 0.5f;

        // source position in continuous space
        const float xs = xd * rx;
        const float ys = yd * ry;

        // integer base index
        const int ig = int(std::floor(xs));
        const int jg = int(std::floor(ys));

        // fractional part
        const float u = xs - float(ig);
        const float v = ys - float(jg);

        // clamp base index so bilinear neighborhood is valid
        const int ig0 = std::clamp(ig, 0, array.shape.x - 2);
        const int jg0 = std::clamp(jg, 0, array.shape.y - 2);

        tile(i, j) = array.get_value_bilinear_at(ig0, jg0, u, v);
      }
  };

  for_each_tile(*this, lambda, cm);
}

void VirtualArray::from_array_bicubic(const Array &array, const ComputeMode &cm)
{
  if (!validate_non_empty(array)) return;

  auto lambda = [&array, this](Array &tile, const TileRegion &region)
  {
    const float rx = float(array.shape.x - 1) / float(this->shape.x - 1);
    const float ry = float(array.shape.y - 1) / float(this->shape.y - 1);

    const glm::ivec2 ij0 = this->tile_region_global_indices(region);

    for (int j = 0; j < region.shape.y; ++j)
      for (int i = 0; i < region.shape.x; ++i)
      {
        // destination pixel center (global)
        const float xd = float(ij0.x + i) + 0.5f;
        const float yd = float(ij0.y + j) + 0.5f;

        // source position in continuous space
        const float xs = xd * rx;
        const float ys = yd * ry;

        // integer base index
        const int ig = int(std::floor(xs));
        const int jg = int(std::floor(ys));

        // fractional part
        const float u = xs - float(ig);
        const float v = ys - float(jg);

        // clamp base index so bilinear neighborhood is valid
        const int ig0 = std::clamp(ig, 0, array.shape.x - 1);
        const int jg0 = std::clamp(jg, 0, array.shape.y - 1);

        tile(i, j) = array.get_value_bicubic_at(ig0, jg0, u, v);
      }
  };

  for_each_tile(*this, lambda, cm);
}

// Access individual cells (slower)
float VirtualArray::get(int global_i, int global_j) const
{
  TileRegion   region = tile_region_from_global_index(global_i, global_j);
  const Array &tile = this->storage.get()->get_tile(region);
  glm::ivec2   local = local_indices(region, global_i, global_j);

  float value = tile(local);
  this->storage->release_tile(region);
  return value;
}

glm::ivec2 VirtualArray::get_max_tiles() const
{
  int nx = ceil_div(this->shape.x, this->tile_shape.x);
  int ny = ceil_div(this->shape.y, this->tile_shape.y);

  return {nx, ny};
}

float VirtualArray::get_bilinear(float x, float y) const
{
  float xr = (x - this->bbox.x) / (this->bbox.y - this->bbox.x);
  float yr = (y - this->bbox.z) / (this->bbox.w - this->bbox.z);

  xr = std::clamp(xr, 0.f, 1.f);
  yr = std::clamp(yr, 0.f, 1.f);

  float xi = xr * (this->shape.x - 1.f);
  float yi = yr * (this->shape.y - 1.f);

  int global_i = int(xi);
  int global_j = int(yi);

  float u = xi - float(global_i);
  float v = yi - float(global_j);

  if (global_i >= this->shape.x - 1)
  {
    global_i = this->shape.x - 2;
    u = 1.f;
  }
  if (global_j >= this->shape.y - 1)
  {
    global_j = this->shape.y - 2;
    v = 1.f;
  }

  int global_i1 = global_i + 1;
  int global_j1 = global_j + 1;

  float v00 = this->get(global_i, global_j);
  float v10 = this->get(global_i1, global_j);
  float v01 = this->get(global_i, global_j1);
  float v11 = this->get(global_i1, global_j1);

  float a10 = v10 - v00;
  float a01 = v01 - v00;
  float a11 = v11 - v10 - v01 + v00;

  return v00 + a10 * u + a01 * v + a11 * u * v;
}

float VirtualArray::get_nearest(float x, float y) const
{
  float xr = (x - this->bbox.x) / (this->bbox.y - this->bbox.x);
  float yr = (y - this->bbox.z) / (this->bbox.w - this->bbox.z);

  xr = std::clamp(xr, 0.f, 1.f);
  yr = std::clamp(yr, 0.f, 1.f);

  int global_i = int(xr * (this->shape.x - 1.f));
  int global_j = int(yr * (this->shape.y - 1.f));

  return this->get(global_i, global_j);
}

int VirtualArray::get_ntiles() const
{
  glm::ivec2 m = this->get_max_tiles();
  return m.x * m.y;
}

std::string VirtualArray::info_string(const ComputeMode &cm, int indent) const
{
  std::ostringstream oss;
  const std::string  pad = std::string(indent, ' ');

  const glm::ivec2 grid = this->get_max_tiles();

  oss << pad << "VirtualArray @" << static_cast<const void *>(this) << "\n";
  oss << pad << "  shape        : " << shape.x << " x " << shape.y << "\n";
  oss << pad << "  tile_size    : " << tile_shape.x << " x " << tile_shape.y
      << "\n";
  oss << pad << "  tile_grid    : " << grid.x << " x " << grid.y << " ("
      << grid.x * grid.y << " tiles)\n";

  oss << pad << "  halo         : " << halo << "\n";

  oss << pad
      << "  storage      : " << (storage ? storage->info_string() : "<none>")
      << "\n";

  oss << pad << "  value_range  : [" << this->min(cm) << ", " << this->max(cm)
      << "]\n";

  return oss.str();
}

glm::ivec2 VirtualArray::local_indices(const TileRegion &region,
                                       int               global_i,
                                       int               global_j) const
{
  // tile indices
  int tx = int(global_i / this->tile_shape.x);
  int ty = int(global_j / this->tile_shape.y);

  // shift from global to local
  int di = tx * this->tile_shape.x - region.halo.x;
  int dj = ty * this->tile_shape.y - region.halo.z;

  glm::ivec2 local = {global_i - di, global_j - dj};

  return local;
}

void VirtualArray::set(int global_i, int global_j, float v)
{
  TileRegion region = tile_region_from_global_index(global_i, global_j);
  Array     &tile = this->storage->get_tile(region);
  glm::ivec2 local = local_indices(region, global_i, global_j);
  tile(local) = v;
  this->storage->release_tile(region);
}

TileRegion VirtualArray::tile_region_from_global_index(int global_i,
                                                       int global_j) const
{
  // belonging tile indices
  int tx = int(global_i / this->tile_shape.x);
  int ty = int(global_j / this->tile_shape.y);

  return this->tile_region_from_tile_coords(tx, ty);
}

TileRegion VirtualArray::tile_region_from_tile_coords(int tile_x,
                                                      int tile_y) const
{
  // actual tile shape
  int w = std::min(this->tile_shape.x,
                   this->shape.x - tile_x * this->tile_shape.x);
  int h = std::min(this->tile_shape.y,
                   this->shape.y - tile_y * this->tile_shape.y);

  // add halo
  glm::ivec4 halo4;

  int tx_max = ceil_div(this->shape.x, this->tile_shape.x) - 1;
  int ty_max = ceil_div(this->shape.y, this->tile_shape.y) - 1;

  halo4 = {(tile_x == 0) ? 0 : this->halo,
           (tile_x == tx_max) ? 0 : this->halo,
           (tile_y == 0) ? 0 : this->halo,
           (tile_y == ty_max) ? 0 : this->halo};

  float fx1 = this->bbox.x + (this->bbox.y - this->bbox.x) *
                                 (tile_x * this->tile_shape.x - halo4.x) /
                                 float(this->shape.x);
  float fx2 = this->bbox.x + (this->bbox.y - this->bbox.x) *
                                 (tile_x * this->tile_shape.x + w + halo4.y) /
                                 float(this->shape.x);
  float fy1 = this->bbox.z + (this->bbox.w - this->bbox.z) *
                                 (tile_y * this->tile_shape.y - halo4.z) /
                                 float(this->shape.y);
  float fy2 = this->bbox.z + (this->bbox.w - this->bbox.z) *
                                 (tile_y * this->tile_shape.y + h + halo4.w) /
                                 float(this->shape.y);

  TileKey key = {.tx = tile_x, .ty = tile_y};

  return TileRegion(key,
                    glm::vec4(fx1, fx2, fy1, fy2),
                    glm::ivec2(w + halo4.x + halo4.y, h + halo4.z + halo4.w),
                    halo4);
}

glm::vec2 VirtualArray::tile_region_global_position(
    const TileRegion &region) const
{
  // tile relative position
  glm::vec2 size = {this->bbox.y - this->bbox.x, this->bbox.w - this->bbox.z};
  glm::vec2 pos = {(region.bbox.x - this->bbox.x) / size.x,
                   (region.bbox.z - this->bbox.z) / size.y};
  return pos;
}

glm::ivec2 VirtualArray::tile_region_global_indices(
    const TileRegion &region) const
{
  int di = region.key.tx * this->tile_shape.x - region.halo.x;
  int dj = region.key.ty * this->tile_shape.y - region.halo.z;
  return {di, dj};
}

Array VirtualArray::to_array(const glm::ivec2  &array_shape,
                             const ComputeMode &cm) const
{
  if (!validate_shape(array_shape)) return Array();

  Array array(array_shape);

  auto lambda = [&array, this](const Array &tile, const TileRegion &region)
  {
    const float rx = float(array.shape.x - 1) / float(this->shape.x - 1);
    const float ry = float(array.shape.y - 1) / float(this->shape.y - 1);

    const glm::ivec2  ij0 = this->tile_region_global_indices(region);
    const glm::ivec4 &b = region.halo;

    // iterate only inner tile cells
    for (int j = b.z; j < region.shape.y - b.w; ++j)
      for (int i = b.x; i < region.shape.x - b.y; ++i)
      {
        // source pixel center (global VA space)
        const float xs = float(ij0.x + i) + 0.5f;
        const float ys = float(ij0.y + j) + 0.5f;

        // destination pixel index
        int ig = int(std::floor(xs * rx));
        int jg = int(std::floor(ys * ry));

        // clamp
        ig = std::clamp(ig, 0, array.shape.x - 1);
        jg = std::clamp(jg, 0, array.shape.y - 1);

        array(ig, jg) = tile(i, j);
      }
  };

  for_each_tile(*this, lambda, cm);
  return array;
}

Array VirtualArray::to_array(const ComputeMode &cm) const
{
  return this->to_array(this->shape, cm);
}

Array VirtualArray::to_array_dbg() const
{
  Array array(this->shape);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      array(i, j) = this->get(i, j);

  return array;
}

void VirtualArray::trim_storage() const
{
  if (this->storage) this->storage->trim();
}

size_t VirtualArray::live_tile_count() const
{
  return this->storage ? this->storage->live_tile_count() : 0;
}

size_t VirtualArray::live_memory_bytes() const
{
  return this->storage ? this->storage->live_memory_bytes() : 0;
}

// --- FUNCTIONS

void copy_data(const VirtualArray &src,
               VirtualArray       &dst,
               const ComputeMode  &cm)
{
  for_each_tile(
      {&src},
      {&dst},
      [](std::vector<const Array *> p_arrays_in,
         std::vector<Array *>       p_arrays_out,
         const TileRegion &) { *p_arrays_out[0] = *p_arrays_in[0]; },
      cm);
}

} // namespace hmap
