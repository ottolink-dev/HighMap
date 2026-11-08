/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file tile_region.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <string>

#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

#include "highmap/algebra.hpp"

namespace hmap
{

struct TileKey
{
  int tx;
  int ty;

  bool operator==(const TileKey &other) const
  {
    return tx == other.tx && ty == other.ty;
  }

  unsigned int hash() const
  {
    return (static_cast<uint64_t>(tx) << 32) | (static_cast<uint32_t>(ty));
  }
};

struct TileRegion
{
  TileKey    key;          // ID
  glm::vec4  bbox;         // (x1, x2, y1, y2) in arbitrary coordinates
  glm::ivec2 shape;        // number of cells in tile (incl. halo)
  glm::ivec4 halo;         // number of extra cells around overlapping with
                           // other
                           // tiles: (west, east, south, north)
  glm::ivec2 global_shape; // total domain shape (nx, ny)
  glm::ivec2 tile_shape;   // nominal tile shape (excluding halo)
  glm::ivec2 tile_grid;    // number of tiles along (x, y)

  TileRegion(const TileKey    &key = {},
             const glm::vec4  &bbox = {0.f, 0.f, 0.f, 0.f},
             const glm::ivec2 &shape = {0, 0},
             const glm::ivec4 &halo = {0, 0, 0, 0},
             const glm::ivec2 &global_shape = {0, 0},
             const glm::ivec2 &tile_shape = {0, 0},
             const glm::ivec2 &tile_grid = {0, 0});

  static const TileRegion &dummy();

  glm::vec2   cell_center(int i, int j) const;
  glm::vec2   cell_corner(int i, int j) const;
  bool        contains_global(int gx, int gy) const;
  bool        contains_global_core(int gx, int gy) const;
  glm::ivec2  global_offset() const;
  glm::ivec2  global_to_local(int gx, int gy) const;
  std::string info_string(int indent = 0) const;
  glm::vec4   inner_bbox() const;
  bool        is_east_boundary() const;
  bool        is_north_boundary() const;
  bool        is_south_boundary() const;
  bool        is_west_boundary() const;
  std::string key_as_string() const;
  int         linear_index() const;
  glm::ivec2  local_to_global(int lx, int ly) const;
};

} // namespace hmap
