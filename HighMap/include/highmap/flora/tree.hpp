/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file tree.hpp
 * @copyright Copyright (c) 2025 Otto Link.
 */
#pragma once

#include <glm/glm.hpp>

#include "highmap/geometry/point.hpp"

namespace hmap
{

/**
 * @class Tree
 * @brief Represents a single tree instance with position, species mark, and
 * radius.
 */
class Tree
{
public:
  glm::vec3 position = {0.f, 0.f, 0.f}; ///< Tree 3D position (x, y, elevation
                                        // z).
  uint32_t species_id = 0;              ///< Species identifier or category
                                        // mark.
  float radius = 1.0f;                  ///< Canopy / collision radius (global
                                        // scale).

  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Default constructor initializing tree at origin.
   */
  Tree() = default;

  /**
   * @brief Constructs a Tree from 3D position, species, and radius.
   * @param position   3D coordinates (x, y, z).
   * @param species_id Species identifier.
   * @param radius     Canopy / collision radius (global scale).
   */
  Tree(const glm::vec3 &position, uint32_t species_id = 0, float radius = 1.0f);

  /**
   * @brief Constructs a Tree from 2D position (z = 0), species, and radius.
   * @param position_2d 2D coordinates (x, y).
   * @param species_id  Species identifier.
   * @param radius      Canopy / collision radius (global scale).
   */
  Tree(const glm::vec2 &position_2d,
       uint32_t         species_id = 0,
       float            radius = 1.0f);

  /**
   * @brief Constructs a Tree from coordinate components.
   * @param x          The x-coordinate.
   * @param y          The y-coordinate.
   * @param z          The elevation / z-coordinate.
   * @param species_id Species identifier.
   * @param radius     Canopy / collision radius (global scale).
   */
  Tree(float    x,
       float    y,
       float    z = 0.0f,
       uint32_t species_id = 0,
       float    radius = 1.0f);

  // ==========================================================================
  //  Conversions
  // ==========================================================================

  /**
   * @brief Converts the tree to a 2D Point (x, y) with value set to radius.
   * @return Point Point representation.
   */
  Point to_point() const;

  /**
   * @brief Converts the tree position to glm::vec2 (x, y).
   * @return glm::vec2 The 2D coordinates.
   */
  glm::vec2 to_vec2() const;

  /**
   * @brief Converts the tree position to glm::vec3 (x, y, z).
   * @return glm::vec3 The 3D coordinates.
   */
  glm::vec3 to_vec3() const;

  // ==========================================================================
  //  Operators
  // ==========================================================================

  bool operator==(const Tree &other) const;
  bool operator!=(const Tree &other) const;
};

} // namespace hmap
