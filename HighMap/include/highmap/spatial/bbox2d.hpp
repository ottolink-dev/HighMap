/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file bbox2d.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Axis-aligned 2D bounding box in world coordinates.
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <algorithm>

#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

namespace hmap
{

/**
 * @brief Axis-aligned 2D bounding box in world coordinates.
 */
struct BBox2D
{
  glm::vec2 min{0.f, 0.f};
  glm::vec2 max{1.f, 1.f};

  constexpr BBox2D() = default;
  constexpr BBox2D(const glm::vec2 &min_p, const glm::vec2 &max_p)
      : min(min_p), max(max_p)
  {
  }
  constexpr BBox2D(float xmin, float xmax, float ymin, float ymax)
      : min(xmin, ymin), max(xmax, ymax)
  {
  }

  static constexpr BBox2D unit()
  {
    return BBox2D(glm::vec2(0.f, 0.f), glm::vec2(1.f, 1.f));
  }

  constexpr bool contains(const glm::vec2 &p) const
  {
    return p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y;
  }

  constexpr glm::vec2 clamp(const glm::vec2 &p) const
  {
    return {std::clamp(p.x, min.x, max.x), std::clamp(p.y, min.y, max.y)};
  }

  constexpr glm::vec2 size() const
  {
    return max - min;
  }

  constexpr glm::vec2 center() const
  {
    return 0.5f * (min + max);
  }

  constexpr glm::vec4 to_vec4() const
  {
    return {min.x, max.x, min.y, max.y};
  }

  static constexpr BBox2D from_vec4(const glm::vec4 &v)
  {
    return {{v.x, v.z}, {v.y, v.w}};
  }
};

} // namespace hmap
