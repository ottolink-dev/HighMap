#include <cstdint>
#include <locale>
#include <sstream>
#include <string>

#include "highmap/geometry/point.hpp"
#include "highmap/scatter/scatter_item.hpp"

namespace hmap
{

// ==========================================================================
//  Constructors
// ==========================================================================

ScatterItem::ScatterItem(const glm::vec3 &position,
                         uint32_t         class_id,
                         float            radius)
    : position(position), class_id(class_id), radius(radius)
{
}

ScatterItem::ScatterItem(const glm::vec2 &position_2d,
                         uint32_t         class_id,
                         float            radius)
    : position(glm::vec3(position_2d.x, position_2d.y, 0.0f)),
      class_id(class_id),
      radius(radius)
{
}

ScatterItem::ScatterItem(float    x,
                         float    y,
                         float    z,
                         uint32_t class_id,
                         float    radius)
    : position(glm::vec3(x, y, z)), class_id(class_id), radius(radius)
{
}

// ==========================================================================
//  Conversions
// ==========================================================================

Point ScatterItem::to_point() const
{
  return Point(position.x, position.y, radius);
}

std::string ScatterItem::to_string() const
{
  std::ostringstream ss;
  ss.imbue(std::locale("C"));
  ss << "ScatterItem(pos=[" << position.x << ", " << position.y << ", "
     << position.z << "], class_id=" << class_id << ", radius=" << radius
     << ")";
  return ss.str();
}

glm::vec2 ScatterItem::to_vec2() const
{
  return glm::vec2(position.x, position.y);
}

glm::vec3 ScatterItem::to_vec3() const
{
  return position;
}

// ==========================================================================
//  Operators
// ==========================================================================

bool ScatterItem::operator==(const ScatterItem &other) const
{
  return position == other.position && class_id == other.class_id &&
         radius == other.radius;
}

bool ScatterItem::operator!=(const ScatterItem &other) const
{
  return !(*this == other);
}

} // namespace hmap
