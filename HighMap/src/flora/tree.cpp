/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include "highmap/flora/tree.hpp"

namespace hmap
{

// ==========================================================================
//  Constructors
// ==========================================================================

Tree::Tree(const glm::vec3 &position, uint32_t species_id, float radius)
    : position(position), species_id(species_id), radius(radius)
{
}

Tree::Tree(const glm::vec2 &position_2d, uint32_t species_id, float radius)
    : position(glm::vec3(position_2d.x, position_2d.y, 0.0f)),
      species_id(species_id),
      radius(radius)
{
}

Tree::Tree(float x, float y, float z, uint32_t species_id, float radius)
    : position(glm::vec3(x, y, z)), species_id(species_id), radius(radius)
{
}

// ==========================================================================
//  Conversions
// ==========================================================================

Point Tree::to_point() const
{
  return Point(position.x, position.y, radius);
}

glm::vec2 Tree::to_vec2() const
{
  return glm::vec2(position.x, position.y);
}

glm::vec3 Tree::to_vec3() const
{
  return position;
}

// ==========================================================================
//  Operators
// ==========================================================================

bool Tree::operator==(const Tree &other) const
{
  return position == other.position && species_id == other.species_id &&
         radius == other.radius;
}

bool Tree::operator!=(const Tree &other) const
{
  return !(*this == other);
}

} // namespace hmap
