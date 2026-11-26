/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file uv_rect.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Normalized domain/UV rectangle [u_min, u_max] x [v_min, v_max].
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <algorithm>

#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

#include "highmap/spatial/bbox2d.hpp"

namespace hmap
{

/**
 * @brief Axis-aligned rectangle in normalized [0, 1]^2 domain/UV space.
 *
 * Represents the parametric sub-region of a tile or Array within the unit
 * terrain domain. (Halo cells may extend slightly outside [0, 1].)
 */
struct UvRect
{
  float u_min = 0.f;
  float u_max = 1.f;
  float v_min = 0.f;
  float v_max = 1.f;

  constexpr UvRect() = default;
  constexpr UvRect(float u0, float u1, float v0, float v1)
      : u_min(u0), u_max(u1), v_min(v0), v_max(v1)
  {
  }
  constexpr UvRect(const glm::vec4 &v)
      : u_min(v.x), u_max(v.y), v_min(v.z), v_max(v.w)
  {
  }

  static constexpr UvRect unit()
  {
    return {0.f, 1.f, 0.f, 1.f};
  }

  static constexpr UvRect from_vec4(const glm::vec4 &v)
  {
    return {v.x, v.y, v.z, v.w};
  }

  constexpr glm::vec4 to_vec4() const
  {
    return {u_min, u_max, v_min, v_max};
  }

  constexpr BBox2D to_bbox2d() const
  {
    return BBox2D(glm::vec2(u_min, v_min), glm::vec2(u_max, v_max));
  }

  static constexpr UvRect from_bbox2d(const BBox2D &bbox)
  {
    return UvRect(bbox.min.x, bbox.max.x, bbox.min.y, bbox.max.y);
  }

  constexpr glm::vec2 size() const
  {
    return {u_max - u_min, v_max - v_min};
  }

  constexpr bool contains(float u, float v) const
  {
    return u >= u_min && u <= u_max && v >= v_min && v <= v_max;
  }

  constexpr glm::vec2 lerp(float rx, float ry) const
  {
    return {u_min + rx * (u_max - u_min), v_min + ry * (v_max - v_min)};
  }

  constexpr glm::vec2 normalize(float u, float v) const
  {
    const float du = u_max - u_min;
    const float dv = v_max - v_min;
    return {(du != 0.f) ? (u - u_min) / du : 0.f,
            (dv != 0.f) ? (v - v_min) / dv : 0.f};
  }
};

} // namespace hmap
