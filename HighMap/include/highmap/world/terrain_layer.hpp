/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file terrain_layer.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief World-space terrain layer pairing a UV-domain storage (Array or
 * VirtualArray) with a CoordFrame and vertical elevation scale/offset.
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <glm/vec2.hpp>

#include "highmap/array.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/spatial/bbox2d.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

/**
 * @brief Pairs a normalized [0, 1]^2 UV-domain height storage (Array or
 * VirtualArray) with a world-space CoordFrame and vertical scale/offset.
 */
template <typename Storage> struct TerrainLayer
{
  Storage    height;
  CoordFrame frame;
  float      z_scale = 1.f;
  float      z_offset = 0.f;

  BBox2D world_bbox() const
  {
    return BBox2D::from_vec4(frame.compute_bounding_box());
  }

  glm::vec2 world_cell_size() const
  {
    const glm::vec2 sz = frame.get_size();
    return {
        (height.shape.x > 0) ? sz.x / static_cast<float>(height.shape.x) : 0.f,
        (height.shape.y > 0) ? sz.y / static_cast<float>(height.shape.y) : 0.f};
  }
};

using TerrainPatch = TerrainLayer<Array>;
using VirtualTerrain = TerrainLayer<VirtualArray>;

} // namespace hmap
