/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file virtual_array.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <atomic>
#include <future>
#include <thread>

#include "macrologger.h"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/math/core.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

namespace hmap
{

// =====================================
// Peripheric classes
// =====================================

enum ForEachMode : int
{
  VA_SEQUENTIAL,             // tile-by-tile, single thread
  VA_DISTRIBUTED,            // tile-by-tile, multi-threaded
  VA_SINGLE_ARRAY,           // full array materialized at once
  VA_SINGLE_ARRAY_DOWNSCALED // full array materialized at once
};

inline std::string to_string(ForEachMode m)
{
  switch (m)
  {
  case ForEachMode::VA_SEQUENTIAL: return "VA_SEQUENTIAL";
  case ForEachMode::VA_DISTRIBUTED: return "VA_DISTRIBUTED";
  case ForEachMode::VA_SINGLE_ARRAY: return "VA_SINGLE_ARRAY";
  case ForEachMode::VA_SINGLE_ARRAY_DOWNSCALED:
    return "VA_SINGLE_ARRAY_DOWNSCALED";
  }
  return "UNKNOWN";
}

static std::map<std::string, int> for_each_mode_as_string = {
    {"Distributed", ForEachMode::VA_DISTRIBUTED},
    {"Sequential", ForEachMode::VA_SEQUENTIAL},
    {"Single array", ForEachMode::VA_SINGLE_ARRAY},
    {"Single array downscaled", ForEachMode::VA_SINGLE_ARRAY_DOWNSCALED}};

struct ComputeMode
{
  ForEachMode mode;
  bool        trim_storage = false;
  int         stride = 1;
  float       k_cutoff = 1.f; // freq. cut-off ratio in ]0, 1]
};

// =====================================
// VirtualArray class
// =====================================

struct VirtualArray
{
  // --- Ctor...

  VirtualArray() = default;

  VirtualArray(glm::ivec2                   shape,
               glm::vec4                    bbox, // (x1, x2, y1, y2)
               glm::ivec2                   tile_shape,
               int                          halo,
               std::unique_ptr<TileStorage> storage);

  VirtualArray(glm::ivec2  shape,
               glm::vec4   bbox,
               glm::ivec2  tile_shape,
               int         halo,
               StorageMode storage_mode = StorageMode::VA_RAM);

  VirtualArray(glm::ivec2  shape,
               glm::ivec2  tile_shape,
               int         halo,
               StorageMode storage_mode = StorageMode::VA_RAM);

  // --- Duplicate

  std::unique_ptr<VirtualArray> clone(const ComputeMode &cm,
                                      bool               deep_copy = false);
  void                          copy_from(VirtualArray      &src,
                                          const ComputeMode &cm,
                                          bool               copy_src_data = true);

  // --- Access individual cells (slower)

  float       get(int global_i, int global_j) const;
  float       get_bilinear(float x, float y) const;
  float       get_nearest(float x, float y) const;
  void        set(int global_i, int global_j, float v);
  glm::ivec2  get_max_tiles() const;
  int         get_ntiles() const;
  std::string info_string(const ComputeMode &cm, int indent = 0) const;

  void  fill(float value, const ComputeMode &cm);
  void  from_array(const Array &array, const ComputeMode &cm);
  void  from_array_bilinear(const Array &array, const ComputeMode &cm);
  void  from_array_bicubic(const Array &array, const ComputeMode &cm);
  Array to_array(const glm::ivec2 &array_shape, const ComputeMode &cm) const;
  Array to_array(const ComputeMode &cm) const;
  Array to_array_dbg() const;

  // --- Processing methods

  float max(const ComputeMode &cm) const;
  float mean(const ComputeMode &cm) const;
  float min(const ComputeMode &cm) const;
  float sum(const ComputeMode &cm) const;

  glm::vec2 range(const ComputeMode &cm) const;
  glm::vec2 range_percentile(float              p_low,
                             float              p_high,
                             const ComputeMode &cm,
                             size_t             bins = 1024) const;

  void inverse(const ComputeMode &cm);
  void remap(float vmin, float vmax, const ComputeMode &cm);
  void remap(float              vmin,
             float              vmax,
             float              from_min,
             float              from_max,
             const ComputeMode &cm);

  std::vector<float> unique_values(const ComputeMode &cm) const;

  void smooth_overlap_buffers();

  // Find which tile covers a given index
  glm::vec2  tile_region_global_position(const TileRegion &region) const;
  glm::ivec2 tile_region_global_indices(const TileRegion &region) const;
  TileRegion tile_region_from_global_index(int global_i, int global_j) const;
  TileRegion tile_region_from_tile_coords(int tile_x, int tile_y) const;
  glm::ivec2 local_indices(const TileRegion &region,
                           int               global_i,
                           int               global_j) const;

  void trim_storage() const;

  // --- Members

  glm::ivec2                   shape;
  glm::vec4                    bbox; // (x1, x2, y1, y2)
  glm::ivec2                   tile_shape;
  int                          halo;
  std::unique_ptr<TileStorage> storage;
};

// actual implementation
#include "highmap/virtual_array/virtual_array.inl"

// functions
void copy_data(VirtualArray &src, VirtualArray &dst, const ComputeMode &cm);

} // namespace hmap
