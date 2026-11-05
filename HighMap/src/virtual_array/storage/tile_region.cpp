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
                       const glm::ivec4 &halo)
    : key(key), bbox(bbox), shape(shape), halo(halo)
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

std::string TileRegion::key_as_string() const
{
  return std::to_string(this->key.tx) + "_" + std::to_string(this->key.ty);
}

} // namespace hmap
