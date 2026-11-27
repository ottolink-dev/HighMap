/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <sstream>
#include <string>

#include <glm/glm.hpp>

#include "highmap/spatial/uv_tile.hpp"

namespace hmap
{

UvTile::UvTile(const TileKey    &key,
               const UvRect     &uv_rect,
               const glm::ivec2 &shape,
               const glm::ivec4 &halo,
               const glm::ivec2 &global_shape,
               const glm::ivec2 &tile_shape,
               const glm::ivec2 &tile_grid)
    : key(key),
      uv_rect(uv_rect),
      shape(shape),
      halo(halo),
      global_shape(global_shape),
      tile_shape(tile_shape),
      tile_grid(tile_grid)
{
}

const UvTile &UvTile::unit()
{
  static const UvTile s_unit;
  return s_unit;
}

const UvTile &UvTile::dummy()
{
  return UvTile::unit();
}

glm::vec2 UvTile::cell_center(int i, int j) const
{
  const float dx = (uv_rect.u_max - uv_rect.u_min) / float(shape.x);
  const float dy = (uv_rect.v_max - uv_rect.v_min) / float(shape.y);

  float x = uv_rect.u_min + (i + 0.5f) * dx;
  float y = uv_rect.v_min + (j + 0.5f) * dy;

  return {x, y};
}

glm::vec2 UvTile::cell_corner(int i, int j) const
{
  float x = glm::mix(uv_rect.u_min, uv_rect.u_max, float(i) / shape.x);
  float y = glm::mix(uv_rect.v_min, uv_rect.v_max, float(j) / shape.y);
  return {x, y};
}

bool UvTile::contains_global(int gx, int gy) const
{
  glm::ivec2 local = global_to_local(gx, gy);

  return (local.x >= 0 && local.x < shape.x && local.y >= 0 &&
          local.y < shape.y);
}

bool UvTile::contains_global_core(int gx, int gy) const
{
  glm::ivec2 local = global_to_local(gx, gy);

  return (local.x >= halo.x && local.x < shape.x - halo.y &&
          local.y >= halo.z && local.y < shape.y - halo.w);
}

glm::ivec2 UvTile::global_offset() const
{
  int gx0 = key.tx * tile_shape.x - halo.x;
  int gy0 = key.ty * tile_shape.y - halo.z;

  return {gx0, gy0};
}

glm::ivec2 UvTile::global_to_local(int gx, int gy) const
{
  glm::ivec2 offset = global_offset();

  return {gx - offset.x, gy - offset.y};
}

std::string UvTile::info_string(int indent) const
{
  std::ostringstream oss;
  const std::string  pad(indent, ' ');

  oss << pad << "UvTile @" << static_cast<const void *>(this) << "\n";
  oss << pad << "  key          : (" << key.tx << ", " << key.ty << ")\n";
  oss << pad << "  shape        : " << shape.x << " x " << shape.y << "\n";
  oss << pad << "  uv_rect      : [" << uv_rect.u_min << ", " << uv_rect.u_max
      << "] x [" << uv_rect.v_min << ", " << uv_rect.v_max << "]\n";
  oss << pad << "  halo         : [" << halo.x << ", " << halo.y << ", "
      << halo.z << ", " << halo.w << "]\n";
  oss << pad << "  tile shape   : " << tile_shape.x << " x " << tile_shape.y
      << "\n";
  oss << pad << "  tile grid    : " << tile_grid.x << " x " << tile_grid.y
      << "\n";
  oss << pad << "  global shape : " << global_shape.x << " x " << global_shape.y
      << "\n";
  oss << pad
      << "  cell size    : " << (uv_rect.u_max - uv_rect.u_min) / float(shape.x)
      << " x " << (uv_rect.v_max - uv_rect.v_min) / float(shape.y) << "\n";

  return oss.str();
}

UvRect UvTile::core_uv_rect() const
{
  const float dx = (uv_rect.u_max - uv_rect.u_min) / float(shape.x);
  const float dy = (uv_rect.v_max - uv_rect.v_min) / float(shape.y);

  return {uv_rect.u_min + float(halo.x) * dx,
          uv_rect.u_max - float(halo.y) * dx,
          uv_rect.v_min + float(halo.z) * dy,
          uv_rect.v_max - float(halo.w) * dy};
}

bool UvTile::is_right_border() const
{
  return halo.y == 0;
}

bool UvTile::is_top_border() const
{
  return halo.w == 0;
}

bool UvTile::is_bottom_border() const
{
  return halo.z == 0;
}

bool UvTile::is_left_border() const
{
  return halo.x == 0;
}

std::string UvTile::key_as_string() const
{
  return std::to_string(this->key.tx) + "_" + std::to_string(this->key.ty);
}

int UvTile::linear_index() const
{
  return key.ty * tile_grid.x + key.tx;
}

glm::ivec2 UvTile::local_to_global(int lx, int ly) const
{
  glm::ivec2 offset = global_offset();

  return {lx + offset.x, ly + offset.y};
}

} // namespace hmap
