/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file uv_tile.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Normalized UV-domain tile descriptor for single-source SPMD Array
 * algorithms and VirtualArray dispatch.
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <cstdint>
#include <string>

#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

#include "highmap/algebra.hpp"
#include "highmap/spatial/uv_rect.hpp"

namespace hmap
{

struct TileKey
{
  int tx = 0;
  int ty = 0;

  bool operator==(const TileKey &other) const
  {
    return tx == other.tx && ty == other.ty;
  }

  unsigned int hash() const
  {
    return (static_cast<uint64_t>(tx) << 32) | (static_cast<uint32_t>(ty));
  }
};

/**
 * @brief Describes an Array tile's position within the normalized [0, 1]^2 UV
 * domain and global index grid.
 */
struct UvTile
{
  TileKey    key;          ///< Tile grid coordinates (tx, ty)
  UvRect     uv_rect;      ///< UV sub-rectangle [u_min, u_max] x [v_min, v_max]
  glm::ivec2 shape;        ///< Number of cells in tile (incl. halo)
  glm::ivec4 halo;         ///< Halo cells: (left, right, bottom, top)
  glm::ivec2 global_shape; ///< Total domain shape (nx, ny)
  glm::ivec2 tile_shape;   ///< Nominal tile shape (excluding halo)
  glm::ivec2 tile_grid;    ///< Number of tiles along (x, y)

  UvTile(const TileKey    &key = {},
         const UvRect     &uv_rect = UvRect::unit(),
         const glm::ivec2 &shape = {0, 0},
         const glm::ivec4 &halo = {0, 0, 0, 0},
         const glm::ivec2 &global_shape = {0, 0},
         const glm::ivec2 &tile_shape = {0, 0},
         const glm::ivec2 &tile_grid = {1, 1});

  UvTile(const UvRect &uv)
      : uv_rect(uv),
        shape({0, 0}),
        halo({0, 0, 0, 0}),
        global_shape({0, 0}),
        tile_shape({0, 0}),
        tile_grid({1, 1})
  {
  }

  static const UvTile &unit();
  static const UvTile &dummy(); ///< @deprecated Use UvTile::unit()

  UvRect core_uv_rect() const;

  glm::vec2 cell_uv_center(int i, int j) const
  {
    return cell_center(i, j);
  }
  glm::vec2 cell_uv_corner(int i, int j) const
  {
    return cell_corner(i, j);
  }
  glm::vec2   cell_center(int i, int j) const;
  glm::vec2   cell_corner(int i, int j) const;
  bool        contains_global(int gx, int gy) const;
  bool        contains_global_core(int gx, int gy) const;
  glm::ivec2  global_offset() const;
  glm::ivec2  global_to_local(int gx, int gy) const;
  std::string info_string(int indent = 0) const;

  bool is_valid() const
  {
    return shape.x > 0 && shape.y > 0;
  }

  bool is_empty() const
  {
    return shape.x <= 0 || shape.y <= 0;
  }

  bool        is_bottom_border() const;
  bool        is_left_border() const;
  bool        is_right_border() const;
  bool        is_top_border() const;
  std::string key_as_string() const;
  int         linear_index() const;
  glm::ivec2  local_to_global(int lx, int ly) const;
};

} // namespace hmap
