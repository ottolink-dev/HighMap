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
#include "highmap/profiling.hpp"
#include "highmap/spatial/uv_tile.hpp"
#include "highmap/virtual_array/tile_storage.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <format>

namespace hmap
{

VirtualArray::VirtualArray(glm::ivec2                   shape,
                           glm::ivec2                   tile_shape,
                           int                          halo,
                           std::unique_ptr<TileStorage> storage)
    : shape(shape),
      tile_shape(tile_shape),
      halo(halo),
      storage(std::move(storage))
{
}

VirtualArray::VirtualArray(glm::ivec2  shape,
                           glm::ivec2  tile_shape,
                           int         halo,
                           StorageMode storage_mode)
    : shape(shape), tile_shape(tile_shape), halo(halo)
{
  this->storage = make_storage(this->shape, this->tile_shape, storage_mode);
}

void VirtualArray::copy_from(const VirtualArray &src,
                             const ComputeMode  &cm,
                             bool                copy_src_data)
{
  if (this == &src) return;

  shape = src.shape;
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
      [value](Array &tile, const UvTile &) { tile = value; },
      cm);
}

void VirtualArray::from_array(const Array &array, const ComputeMode &cm)
{
  HMAP_PROFILE_ZONE("HighMap::VirtualArray::from_array");
  if (!validate_non_empty(array)) return;

  auto lambda = [&array, this](Array &tile, const UvTile &region)
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
  HMAP_PROFILE_ZONE("HighMap::VirtualArray::from_array_bilinear");
  if (!validate_non_empty(array)) return;

  auto lambda = [&array, this](Array &tile, const UvTile &region)
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
  HMAP_PROFILE_ZONE("HighMap::VirtualArray::from_array_bicubic");
  if (!validate_non_empty(array)) return;

  auto lambda = [&array, this](Array &tile, const UvTile &region)
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
  UvTile       region = tile_region_from_global_index(global_i, global_j);
  const Array &tile = this->storage.get()->get_tile(region);
  glm::ivec2   local = local_indices(region, global_i, global_j);

  float value = tile(local);
  this->storage->release_tile(region);
  return value;
}

glm::ivec2 VirtualArray::get_max_tiles() const
{
  if (this->tile_shape.x <= 0 || this->tile_shape.y <= 0 ||
      this->shape.x <= 0 || this->shape.y <= 0)
  {
    return {0, 0};
  }

  int nx = ceil_div(this->shape.x, this->tile_shape.x);
  int ny = ceil_div(this->shape.y, this->tile_shape.y);

  return {nx, ny};
}

float VirtualArray::get_bilinear(float x, float y) const
{
  float xr = std::clamp(x, 0.f, 1.f);
  float yr = std::clamp(y, 0.f, 1.f);

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
  float xr = std::clamp(x, 0.f, 1.f);
  float yr = std::clamp(y, 0.f, 1.f);

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

glm::ivec2 VirtualArray::local_indices(const UvTile &region,
                                       int           global_i,
                                       int           global_j) const
{
  return region.global_to_local(global_i, global_j);
}

void VirtualArray::set(int global_i, int global_j, float v)
{
  UvTile     region = tile_region_from_global_index(global_i, global_j);
  Array     &tile = this->storage->get_tile(region);
  glm::ivec2 local = local_indices(region, global_i, global_j);
  tile(local) = v;
  this->storage->release_tile(region);
}

UvTile VirtualArray::tile_region_from_global_index(int global_i,
                                                   int global_j) const
{
  // belonging tile indices
  int tx = int(global_i / this->tile_shape.x);
  int ty = int(global_j / this->tile_shape.y);

  return this->tile_region_from_tile_coords(tx, ty);
}

UvTile VirtualArray::tile_region_from_tile_coords(int tile_x, int tile_y) const
{
  // actual tile shape
  int w = std::min(this->tile_shape.x,
                   this->shape.x - tile_x * this->tile_shape.x);
  int h = std::min(this->tile_shape.y,
                   this->shape.y - tile_y * this->tile_shape.y);

  // add halo
  glm::ivec4 halo4;

  int nx = ceil_div(this->shape.x, this->tile_shape.x);
  int ny = ceil_div(this->shape.y, this->tile_shape.y);

  int tx_max = nx - 1;
  int ty_max = ny - 1;

  halo4 = {(tile_x == 0) ? 0 : this->halo,
           (tile_x == tx_max) ? 0 : this->halo,
           (tile_y == 0) ? 0 : this->halo,
           (tile_y == ty_max) ? 0 : this->halo};

  float fx1 = (tile_x * this->tile_shape.x - halo4.x) / float(this->shape.x);
  float fx2 = (tile_x * this->tile_shape.x + w + halo4.y) /
              float(this->shape.x);
  float fy1 = (tile_y * this->tile_shape.y - halo4.z) / float(this->shape.y);
  float fy2 = (tile_y * this->tile_shape.y + h + halo4.w) /
              float(this->shape.y);

  TileKey key = {.tx = tile_x, .ty = tile_y};

  return UvTile(key,
                UvRect{fx1, fx2, fy1, fy2},
                glm::ivec2(w + halo4.x + halo4.y, h + halo4.z + halo4.w),
                halo4,
                this->shape,
                this->tile_shape,
                glm::ivec2(nx, ny));
}

glm::vec2 VirtualArray::tile_region_global_position(const UvTile &region) const
{
  return {region.uv_rect.u_min, region.uv_rect.v_min};
}

glm::ivec2 VirtualArray::tile_region_global_indices(const UvTile &region) const
{
  return region.global_offset();
}

Array VirtualArray::to_array(const glm::ivec2  &array_shape,
                             const ComputeMode &cm) const
{
  HMAP_PROFILE_ZONE("HighMap::VirtualArray::to_array");
  if (!validate_shape(array_shape)) return Array();

  Array array(array_shape);

  if (array_shape == this->shape)
  {
    auto lambda = [&array, this](const Array &tile, const UvTile &region)
    {
      const glm::ivec2  ij0 = this->tile_region_global_indices(region);
      const glm::ivec4 &b = region.halo;

      // iterate only inner tile cells
      for (int j = b.z; j < region.shape.y - b.w; ++j)
        for (int i = b.x; i < region.shape.x - b.y; ++i)
          array(ij0.x + i, ij0.y + j) = tile(i, j);
    };

    for_each_tile(*this, lambda, cm);
    return array;
  }

  auto lambda = [&array, this](const Array &tile, const UvTile &region)
  {
    const glm::ivec2  ij0 = this->tile_region_global_indices(region);
    const glm::ivec4 &b = region.halo;

    // global inner-tile bounds in source space
    const int x0 = ij0.x + b.x;
    const int x1 = ij0.x + region.shape.x - b.y;
    const int y0 = ij0.y + b.z;
    const int y1 = ij0.y + region.shape.y - b.w;

    const float rx = (array.shape.x > 1)
                         ? float(this->shape.x - 1) / float(array.shape.x - 1)
                         : 0.f;
    const float ry = (array.shape.y > 1)
                         ? float(this->shape.y - 1) / float(array.shape.y - 1)
                         : 0.f;

    // destination index range whose source coordinates fall in [x0, x1) x [y0,
    // y1)
    int ig0 = (rx > 0.f) ? int(std::ceil(float(x0) / rx)) : 0;
    int ig1 = (rx > 0.f) ? int(std::ceil(float(x1) / rx)) : array.shape.x;
    int jg0 = (ry > 0.f) ? int(std::ceil(float(y0) / ry)) : 0;
    int jg1 = (ry > 0.f) ? int(std::ceil(float(y1) / ry)) : array.shape.y;

    ig0 = std::clamp(ig0, 0, array.shape.x);
    ig1 = std::clamp(ig1, 0, array.shape.x);
    jg0 = std::clamp(jg0, 0, array.shape.y);
    jg1 = std::clamp(jg1, 0, array.shape.y);

    // include the last source pixel in the last tile
    if (x1 == this->shape.x) ig1 = array.shape.x;
    if (y1 == this->shape.y) jg1 = array.shape.y;

    for (int jg = jg0; jg < jg1; ++jg)
    {
      const float ys = std::clamp(float(jg) * ry - float(ij0.y),
                                  0.f,
                                  float(region.shape.y - 1));
      const int   j_lo = int(std::floor(ys));
      const int   j_hi = std::min(j_lo + 1, region.shape.y - 1);
      const float ty = ys - float(j_lo);

      for (int ig = ig0; ig < ig1; ++ig)
      {
        const float xs = std::clamp(float(ig) * rx - float(ij0.x),
                                    0.f,
                                    float(region.shape.x - 1));
        const int   i_lo = int(std::floor(xs));
        const int   i_hi = std::min(i_lo + 1, region.shape.x - 1);
        const float tx = xs - float(i_lo);

        const float v00 = tile(i_lo, j_lo);
        const float v10 = tile(i_hi, j_lo);
        const float v01 = tile(i_lo, j_hi);
        const float v11 = tile(i_hi, j_hi);

        const float v0 = v00 + tx * (v10 - v00);
        const float v1 = v01 + tx * (v11 - v01);
        array(ig, jg) = v0 + ty * (v1 - v0);
      }
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
         const UvTile &) { *p_arrays_out[0] = *p_arrays_in[0]; },
      cm);
}

void sum_per_tile(const VirtualArray &va,
                  std::vector<float> &tile_sums,
                  float              *p_total_inner_sum,
                  const VirtualArray *p_exclusion_mask,
                  float               exclusion_threshold,
                  const ComputeMode  &cm)
{
  // --- Check valid input

  glm::ivec2 num_tiles = va.get_max_tiles();
  int        total_tiles = num_tiles.x * num_tiles.y;

  tile_sums.assign(total_tiles, 0.f);

  if (p_total_inner_sum)
  {
    *p_total_inner_sum = 0.f;
  }

  if (va.empty() || total_tiles == 0)
  {
    return;
  }

  // --- Validate exclusion mask if present

  bool has_exclusion = (p_exclusion_mask != nullptr) &&
                       (p_exclusion_mask->storage != nullptr) &&
                       (p_exclusion_mask->shape.x > 0 &&
                        p_exclusion_mask->shape.y > 0);

  if (has_exclusion)
  {
    if (p_exclusion_mask->shape != va.shape ||
        p_exclusion_mask->tile_shape != va.tile_shape)
    {
      hmap::log::warn(
          "Exclusion VirtualArray shape mismatch with VirtualArray");
      return;
    }
  }

  const VirtualArray *p_excl_va = has_exclusion ? p_exclusion_mask : nullptr;

  float total_inner_sum = 0.f;

  // --- Compute sums per tile

  hmap::for_each_tile(
      {&va, p_excl_va},
      {},
      [&](const std::vector<const Array *> &p_arrays_in,
          std::vector<Array *> &,
          const UvTile &region)
      {
        const auto [p_tile, p_excl_tile] = unpack<2>(p_arrays_in);

        int tile_index = region.key.ty * num_tiles.x + region.key.tx;

        float tile_full_sum = 0.f;
        float tile_inner_sum = 0.f;

        const glm::ivec4 &halo = region.halo;
        int               shape_x = region.shape.x;
        int               shape_y = region.shape.y;

        for (int j = 0; j < shape_y; ++j)
        {
          bool is_inner_y = (j >= halo.z && j < shape_y - halo.w);

          for (int i = 0; i < shape_x; ++i)
          {
            float val = (*p_tile)(i, j);

            if (p_excl_tile && (*p_excl_tile)(i, j) > exclusion_threshold)
            {
              val = 0.0f;
            }

            if (val > 0.0f)
            {
              tile_full_sum += val;

              if (is_inner_y && i >= halo.x && i < shape_x - halo.y)
              {
                tile_inner_sum += val;
              }
            }
          }
        }

        tile_sums[tile_index] = tile_full_sum;
        total_inner_sum += tile_inner_sum;
      },
      cm);

  if (p_total_inner_sum)
  {
    *p_total_inner_sum = total_inner_sum;
  }
}

} // namespace hmap
