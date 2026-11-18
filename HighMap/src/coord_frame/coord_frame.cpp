/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "highmap/coord_frame.hpp"

namespace hmap
{

glm::vec2 AffineTransform2D::apply(float x, float y) const
{
  return {this->ax * x + this->bx * y + this->cx,
          this->ay * x + this->by * y + this->cy};
}

CoordFrame::CoordFrame()
    : origin(glm::vec2(0.f, 0.f)), size(glm::vec2(1.f, 1.f))
{
  this->set_rotation_angle(0.f);
}

CoordFrame::CoordFrame(glm::vec2 origin, glm::vec2 size, float rotation_angle)
    : origin(origin), size(size)
{
  this->set_rotation_angle(rotation_angle);
}

glm::vec4 CoordFrame::compute_bounding_box() const
{
  // Define the four corners relative to the origin
  std::pair<float, float> corners[4] = {{0.f, 0.f},
                                        {1.f, 0.f},
                                        {1.f, 1.f},
                                        {0.f, 1.f}};

  // Rotate each corner and find min/max x and y
  float min_x = std::numeric_limits<float>::max();
  float max_x = std::numeric_limits<float>::lowest();
  float min_y = std::numeric_limits<float>::max();
  float max_y = std::numeric_limits<float>::lowest();

  for (const auto &corner : corners)
  {
    glm::vec2 g = this->map_to_global_coords(corner.first, corner.second);

    // Update bounding box extents
    min_x = std::min(min_x, g.x);
    max_x = std::max(max_x, g.x);
    min_y = std::min(min_y, g.y);
    max_y = std::max(max_y, g.y);
  }

  return glm::vec4(min_x, max_x, min_y, max_y);
}

AffineTransform2D CoordFrame::compute_relative_transform(
    const CoordFrame &other) const
{
  // Maps relative coords (rx, ry) in `this` frame -> global -> relative coords
  // in `other` frame:
  //   g = this->origin + R_this * (size_this * r)
  //   r_other = (1 / size_other) * R_other^T * (g - other.origin)
  const float sx = this->size.x;
  const float sy = this->size.y;

  const float inv_ox = (other.size.x != 0.f) ? (1.f / other.size.x) : 0.f;
  const float inv_oy = (other.size.y != 0.f) ? (1.f / other.size.y) : 0.f;

  // R_other^T * R_this = R(this.angle - other.angle)
  const float cos_rel = this->cos_angle * other.cos_angle +
                        this->sin_angle * other.sin_angle;
  const float sin_rel = this->sin_angle * other.cos_angle -
                        this->cos_angle * other.sin_angle;

  const glm::vec2 d = this->origin - other.origin;
  const float     tx = (d.x * other.cos_angle + d.y * other.sin_angle) * inv_ox;
  const float ty = (-d.x * other.sin_angle + d.y * other.cos_angle) * inv_oy;

  AffineTransform2D tf;
  tf.ax = sx * cos_rel * inv_ox;
  tf.bx = -sy * sin_rel * inv_ox;
  tf.cx = tx;

  tf.ay = sx * sin_rel * inv_oy;
  tf.by = sy * cos_rel * inv_oy;
  tf.cy = ty;

  return tf;
}

glm::vec2 CoordFrame::get_origin() const
{
  return this->origin;
}

glm::vec2 CoordFrame::get_size() const
{
  return this->size;
}

float CoordFrame::get_rotation_angle() const
{
  return this->rotation_angle;
}

void CoordFrame::set_origin(glm::vec2 new_origin)
{
  this->origin = new_origin;
}

void CoordFrame::set_size(glm::vec2 new_size)
{
  this->size = new_size;
}

// Method to determine if a point (x, y) is within the rotated terrain
bool CoordFrame::intersects(const CoordFrame &other) const
{
  const glm::vec4 b1 = this->compute_bounding_box();
  const glm::vec4 b2 = other.compute_bounding_box();

  return !(b1.y < b2.x || b2.y < b1.x || b1.w < b2.z || b2.w < b1.z);
}

bool CoordFrame::is_point_within(float gx, float gy) const
{
  glm::vec2 rel = this->map_to_relative_coords(gx, gy);

  // Check if the unrotated point is within the unrotated bounds
  return this->is_unit_point_within(rel.x, rel.y);
}

bool CoordFrame::is_unit_point_within(float rx, float ry) const
{
  return (rx >= 0.f && rx <= 1.f && ry >= 0.f && ry <= 1.f);
}

// Method to set the rotation angle and update cos_angle and sin_angle
void CoordFrame::set_rotation_angle(float new_angle)
{
  this->rotation_angle = new_angle;

  float angle_rad = this->rotation_angle * M_PI / 180.0f;

  this->cos_angle = std::cos(angle_rad);
  this->sin_angle = std::sin(angle_rad);
}

glm::vec2 CoordFrame::map_to_global_coords(float rx, float ry) const
{
  rx *= this->size.x;
  ry *= this->size.y;

  float gx = this->origin.x + rx * this->cos_angle - ry * this->sin_angle;
  float gy = this->origin.y + rx * this->sin_angle + ry * this->cos_angle;

  return glm::vec2(gx, gy);
}

// Method to apply inverse rotation to a point relative to the terrain origin
glm::vec2 CoordFrame::map_to_relative_coords(float gx, float gy) const
{
  // Translate the point to be relative to the terrain's origin
  glm::vec2 translated = glm::vec2(gx, gy) - this->origin;

  // Unrotate the point using precomputed cos_angle and sin_angle
  float unrotated_x = translated.x * this->cos_angle +
                      translated.y * this->sin_angle;
  float unrotated_y = -translated.x * this->sin_angle +
                      translated.y * this->cos_angle;

  unrotated_x /= this->size.x;
  unrotated_y /= this->size.y;

  return glm::vec2(unrotated_x, unrotated_y);
}

glm::vec2 CoordFrame::map_to_relative_coords(const CoordFrame &other,
                                             float             rx,
                                             float             ry) const
{
  return this->compute_relative_transform(other).apply(rx, ry);
}

float CoordFrame::normalized_distance_to_edges(float gx, float gy) const
{
  // switch to unit-square coordinates
  glm::vec2 rel = this->map_to_relative_coords(gx, gy);

  // distances to the 4 edges, times 2 to get something in [0, 1]
  float dmin = 2.f * std::min(1.f - rel.y,
                              std::min(rel.y, std::min(rel.x, 1.f - rel.x)));

  return dmin;
}

float CoordFrame::normalized_shape_factor(float gx, float gy) const
{
  // switch to unit-square coordinates
  glm::vec2 rel = this->map_to_relative_coords(gx, gy);

  return this->normalized_shape_factor_rel(rel.x, rel.y);
}

float CoordFrame::normalized_shape_factor_rel(float rx, float ry) const
{
  return 256.f * rx * rx * (1.f - rx) * (1.f - rx) * ry * ry * (1.f - ry) *
         (1.f - ry);
}

} // namespace hmap
