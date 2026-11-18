/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file coord_frame.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief
 *
 * @copyright Copyright (c) 2023
 */
#pragma once
#include <cmath>
#include <map>
#include <vector>

#include "highmap/algebra.hpp"

namespace hmap
{

/**
 * @brief Lightweight 2D affine transformation mapping (x, y) -> (ax*x + bx*y +
 * cx, ay*x + by*y + cy).
 */
struct AffineTransform2D
{
  float ax = 1.f;
  float bx = 0.f;
  float cx = 0.f;
  float ay = 0.f;
  float by = 1.f;
  float cy = 0.f;

  glm::vec2 apply(float x, float y) const;
};

class CoordFrame
{
public:
  CoordFrame();

  CoordFrame(glm::vec2 origin, glm::vec2 size, float rotation_angle);

  // Getters
  glm::vec2 get_origin() const;

  glm::vec2 get_size() const;

  float get_rotation_angle() const;

  // Setters
  void set_origin(glm::vec2 new_origin);

  void set_size(glm::vec2 new_size);

  void set_rotation_angle(float new_angle);

  // Methods

  glm::vec4 compute_bounding_box() const;

  bool intersects(const CoordFrame &other) const;

  bool is_point_within(float gx, float gy) const;

  bool is_unit_point_within(float rx, float ry) const;

  AffineTransform2D compute_relative_transform(const CoordFrame &other) const;

  glm::vec2 map_to_global_coords(float rx, float ry) const;

  // relative coords always in [0, 1] x [0, 1] (unit square)
  glm::vec2 map_to_relative_coords(float gx, float gy) const;

  glm::vec2 map_to_relative_coords(const CoordFrame &other,
                                   float             rx,
                                   float             ry) const;

  float normalized_distance_to_edges(float gx, float gy) const;

  float normalized_shape_factor(float gx, float gy) const;

  float normalized_shape_factor_rel(float rx, float ry) const;

private:
  glm::vec2 origin;
  glm::vec2 size;
  float     rotation_angle;

  float cos_angle;
  float sin_angle;
};

} // namespace hmap
