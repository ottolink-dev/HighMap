/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <sstream>
#include <string>

#include <glm/glm.hpp>

#include "highmap/virtual_array/tile_region.hpp"

namespace hmap
{

TileRegion::TileRegion(const TileKey    &key,
                       const glm::vec4  &bbox,
                       const glm::ivec2 &shape,
                       const glm::ivec4 &halo,
                       const glm::ivec2 &global_shape,
                       const glm::ivec2 &tile_shape,
                       const glm::ivec2 &tile_grid)
    : key(key),
      bbox(bbox),
      shape(shape),
      halo(halo),
      global_shape(global_shape),
      tile_shape(tile_shape),
      tile_grid(tile_grid)
{
}

const TileRegion &TileRegion::dummy()
{
  static const TileRegion s_dummy;
  return s_dummy;
}

glm::vec2 TileRegion::cell_center(int i, int j) const
{
  const float x_min = bbox.x;
  const float x_max = bbox.y;
  const float y_min = bbox.z;
  const float y_max = bbox.w;

  const float dx = (x_max - x_min) / float(shape.x);
  const float dy = (y_max - y_min) / float(shape.y);

  float x = x_min + (i + 0.5f) * dx;
  float y = y_min + (j + 0.5f) * dy;

  return {x, y};
}

glm::vec2 TileRegion::cell_corner(int i, int j) const
{
  float x = glm::mix(bbox.x, bbox.y, float(i) / shape.x);
  float y = glm::mix(bbox.z, bbox.w, float(j) / shape.y);
  return {x, y};
}

bool TileRegion::contains_global(int gx, int gy) const
{
  glm::ivec2 local = global_to_local(gx, gy);

  return (local.x >= 0 && local.x < shape.x && local.y >= 0 &&
          local.y < shape.y);
}

bool TileRegion::contains_global_core(int gx, int gy) const
{
  glm::ivec2 local = global_to_local(gx, gy);

  return (local.x >= halo.x && local.x < shape.x - halo.y &&
          local.y >= halo.z && local.y < shape.y - halo.w);
}

glm::ivec2 TileRegion::global_offset() const
{
  int gx0 = key.tx * tile_shape.x - halo.x;
  int gy0 = key.ty * tile_shape.y - halo.z;

  return {gx0, gy0};
}

glm::ivec2 TileRegion::global_to_local(int gx, int gy) const
{
  glm::ivec2 offset = global_offset();

  return {gx - offset.x, gy - offset.y};
}

std::string TileRegion::info_string(int indent) const
{
  std::ostringstream oss;
  const std::string  pad(indent, ' ');

  oss << pad << "TileRegion @" << static_cast<const void *>(this) << "\n";
  oss << pad << "  key          : (" << key.tx << ", " << key.ty << ")\n";
  oss << pad << "  shape        : " << shape.x << " x " << shape.y << "\n";
  oss << pad << "  bbox         : [" << bbox.x << ", " << bbox.y << "] x ["
      << bbox.z << ", " << bbox.w << "]\n";
  oss << pad << "  halo         : [" << halo.x << ", " << halo.y << ", "
      << halo.z << ", " << halo.w << "]\n";
  oss << pad << "  tile shape   : " << tile_shape.x << " x " << tile_shape.y
      << "\n";
  oss << pad << "  tile grid    : " << tile_grid.x << " x " << tile_grid.y
      << "\n";
  oss << pad << "  global shape : " << global_shape.x << " x " << global_shape.y
      << "\n";
  oss << pad << "  cell size    : " << (bbox.y - bbox.x) / float(shape.x)
      << " x " << (bbox.w - bbox.z) / float(shape.y) << "\n";

  return oss.str();
}

glm::vec4 TileRegion::inner_bbox() const
{
  const float dx = (bbox.y - bbox.x) / float(shape.x);
  const float dy = (bbox.w - bbox.z) / float(shape.y);

  float inner_x1 = bbox.x + float(halo.x) * dx;
  float inner_x2 = bbox.y - float(halo.y) * dx;
  float inner_y1 = bbox.z + float(halo.z) * dy;
  float inner_y2 = bbox.w - float(halo.w) * dy;

  return {inner_x1, inner_x2, inner_y1, inner_y2};
}

bool TileRegion::is_east_boundary() const
{
  return halo.y == 0;
}

bool TileRegion::is_north_boundary() const
{
  return halo.w == 0;
}

bool TileRegion::is_south_boundary() const
{
  return halo.z == 0;
}

bool TileRegion::is_west_boundary() const
{
  return halo.x == 0;
}

std::string TileRegion::key_as_string() const
{
  return std::to_string(this->key.tx) + "_" + std::to_string(this->key.ty);
}

int TileRegion::linear_index() const
{
  return key.ty * tile_grid.x + key.tx;
}

glm::ivec2 TileRegion::local_to_global(int lx, int ly) const
{
  glm::ivec2 offset = global_offset();

  return {lx + offset.x, ly + offset.y};
}

} // namespace hmap
