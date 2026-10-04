/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file scatter_item.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstdint>
#include <string>

#include <glm/glm.hpp>

#include "highmap/geometry/point.hpp"

#define HMAP_DEFAULT_SCATTER_RADIUS 1e-3f

namespace hmap
{

/**
 * @class ScatterItem
 * @brief Represents a single spatial instance (tree, rock, prop) with position,
 * class identifier, and radius.
 */
class ScatterItem
{
public:
  glm::vec3 position = {0.f,
                        0.f,
                        0.f}; ///< 3D spatial position (x, y, elevation z).
  uint32_t  class_id =
      0; ///< Unique category/class identifier (species, rock type, prop class).
  float radius = HMAP_DEFAULT_SCATTER_RADIUS; ///< Characteristic radius.

  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Default constructor initializing item at origin.
   */
  ScatterItem() = default;

  /**
   * @brief Constructs a ScatterItem from 3D position, class_id, and radius.
   * @param position 3D coordinates (x, y, z).
   * @param class_id Class / category identifier.
   * @param radius   Characteristic radius.
   */
  ScatterItem(const glm::vec3 &position,
              uint32_t         class_id = 0,
              float            radius = HMAP_DEFAULT_SCATTER_RADIUS);

  /**
   * @brief Constructs a ScatterItem from 2D position (z = 0), class_id, and
   * radius.
   * @param position_2d 2D coordinates (x, y).
   * @param class_id    Class / category identifier.
   * @param radius      Characteristic radius.
   */
  ScatterItem(const glm::vec2 &position_2d,
              uint32_t         class_id = 0,
              float            radius = HMAP_DEFAULT_SCATTER_RADIUS);

  /**
   * @brief Constructs a ScatterItem from coordinate components.
   * @param x        The x-coordinate.
   * @param y        The y-coordinate.
   * @param z        The elevation / z-coordinate.
   * @param class_id Class / category identifier.
   * @param radius   Characteristic radius.
   */
  ScatterItem(float    x,
              float    y,
              float    z = 0.0f,
              uint32_t class_id = 0,
              float    radius = HMAP_DEFAULT_SCATTER_RADIUS);

  // ==========================================================================
  //  Conversions
  // ==========================================================================

  /**
   * @brief Converts the item to a 2D Point (x, y) with value set to radius.
   * @return Point Point representation.
   */
  Point to_point() const;

  /**
   * @brief Returns a formatted string representation of the item.
   * @return std::string Formatted string with item properties.
   */
  std::string to_string() const;

  /**
   * @brief Converts the item position to glm::vec2 (x, y).
   * @return glm::vec2 The 2D coordinates.
   */
  glm::vec2 to_vec2() const;

  /**
   * @brief Converts the item position to glm::vec3 (x, y, z).
   * @return glm::vec3 The 3D coordinates.
   */
  glm::vec3 to_vec3() const;

  // ==========================================================================
  //  Operators
  // ==========================================================================

  bool operator==(const ScatterItem &other) const;
  bool operator!=(const ScatterItem &other) const;
};

} // namespace hmap
