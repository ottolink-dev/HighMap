/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file point.hpp
 * @copyright Copyright (c) 2023 Otto Link.
 */
#pragma once
#include <cmath>
#include <optional>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @class Point
 * @brief A class to represent and manipulate 2D points that can carry a value.
 *
 * See unit tests: @ref test_point.cpp
 */
class Point
{
public:
  float x; ///< The x-coordinate of the point.
  float y; ///< The y-coordinate of the point.
  float v; ///< The value at the point.

  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Default constructor initializing the point to (0, 0, 0).
   */
  constexpr Point() noexcept : x(0.f), y(0.f), v(0.f)
  {
  }

  /**
   * @brief Parameterized constructor initializing the point to given values.
   * @param x The x-coordinate of the point.
   * @param y The y-coordinate of the point.
   * @param v The value at the point.
   */
  constexpr Point(float x, float y, float v = 0.f) noexcept : x(x), y(y), v(v)
  {
  }

  /**
   * @brief Constructs a Point from a glm::vec2 and an optional value.
   * @param xy 2D coordinate vector.
   * @param v  The value at the point.
   */
  constexpr explicit Point(const glm::vec2 &xy, float v = 0.f) noexcept
      : x(xy.x), y(xy.y), v(v)
  {
  }

  /**
   * @brief Constructs a Point from a glm::vec3 (x, y, v).
   * @param xyv 3D vector containing (x, y, v).
   */
  constexpr explicit Point(const glm::vec3 &xyv) noexcept
      : x(xyv.x), y(xyv.y), v(xyv.z)
  {
  }

  // ==========================================================================
  //  Conversions
  // ==========================================================================

  /**
   * @brief Converts the point's coordinates to glm::vec2.
   * @return glm::vec2 The (x, y) coordinates.
   */
  constexpr glm::vec2 to_vec2() const noexcept
  {
    return glm::vec2(x, y);
  }

  /**
   * @brief Converts the point's coordinates and value to glm::vec3.
   * @return glm::vec3 The (x, y, v) components.
   */
  constexpr glm::vec3 to_vec3() const noexcept
  {
    return glm::vec3(x, y, v);
  }

  // ==========================================================================
  //  Operators
  // ==========================================================================

  /**
   * @brief Equality operator to check if two points are the same.
   *
   * @param  other The point to compare with.
   * @return       true if the points are equal, false otherwise.
   */
  constexpr bool operator==(const Point &other) const noexcept
  {
    return (x == other.x && y == other.y && v == other.v);
  }

  /**
   * @brief Inequality operator to check if two points are different.
   *
   * @param  other The point to compare with.
   * @return       true if the points are not equal, false otherwise.
   */
  constexpr bool operator!=(const Point &other) const noexcept
  {
    return !(*this == other);
  }

  /**
   * @brief Adds two points.
   * @param  other The point to add.
   * @return       The result of adding the two points.
   */
  constexpr Point operator+(const Point &other) const noexcept
  {
    return Point(x + other.x, y + other.y, v + other.v);
  }

  /**
   * @brief Subtracts two points.
   * @param  other The point to subtract.
   * @return       The result of subtracting the other point from this point.
   */
  constexpr Point operator-(const Point &other) const noexcept
  {
    return Point(x - other.x, y - other.y, v - other.v);
  }

  /**
   * @brief Multiplies the point by a scalar.
   * @param  scalar The scalar to multiply by.
   * @return        The result of the multiplication.
   */
  constexpr Point operator*(float scalar) const noexcept
  {
    return Point(x * scalar, y * scalar, v * scalar);
  }

  /**
   * @brief Divides the point by a scalar.
   * @param  scalar The scalar to divide by.
   * @return        The result of the division.
   */
  constexpr Point operator/(float scalar) const noexcept
  {
    return Point(x / scalar, y / scalar, v / scalar);
  }

  /**
   * @brief Adds another point to this point in-place.
   * @param  other The point to add.
   * @return       Reference to this point.
   */
  constexpr Point &operator+=(const Point &other) noexcept
  {
    x += other.x;
    y += other.y;
    v += other.v;
    return *this;
  }

  /**
   * @brief Subtracts another point from this point in-place.
   * @param  other The point to subtract.
   * @return       Reference to this point.
   */
  constexpr Point &operator-=(const Point &other) noexcept
  {
    x -= other.x;
    y -= other.y;
    v -= other.v;
    return *this;
  }

  /**
   * @brief Multiplies this point by a scalar in-place.
   * @param  scalar The scalar multiplier.
   * @return        Reference to this point.
   */
  constexpr Point &operator*=(float scalar) noexcept
  {
    x *= scalar;
    y *= scalar;
    v *= scalar;
    return *this;
  }

  /**
   * @brief Divides this point by a scalar in-place.
   * @param  scalar The scalar divisor.
   * @return        Reference to this point.
   */
  constexpr Point &operator/=(float scalar) noexcept
  {
    x /= scalar;
    y /= scalar;
    v /= scalar;
    return *this;
  }

  /**
   * @brief Scalar multiplication (scalar * Point).
   *
   * @param  scalar The scalar value to multiply with.
   * @param  point  The point to multiply.
   * @return        Point A new point with each component multiplied by scalar.
   */
  friend constexpr Point operator*(float scalar, const Point &point) noexcept
  {
    return Point(scalar * point.x, scalar * point.y, scalar * point.v);
  }

  // ==========================================================================
  //  Methods
  // ==========================================================================

  /**
   * @brief Prints the coordinates and value of the Point object.
   *
   * This function outputs the Point's x, y coordinates, and value `v` to stdout
   * in the format `(x, y, v)`, followed by a newline.
   */
  void print() const;

  /**
   * @brief Updates the point's value based on bilinear interpolation from an
   * array.
   *
   * @param array The input `Array` from which the value is interpolated.
   * @param bbox  Bounding box used for normalizing coordinates {xmin, xmax,
   * ymin, ymax}.
   */
  void set_value_from_array(const Array &array, const glm::vec4 &bbox);
};

// ==========================================================================
//  Free Functions
// ==========================================================================

/**
 * @brief Computes the angle between two points relative to the x-axis.
 *
 * @param  p1 The starting point of the vector.
 * @param  p2 The ending point of the vector.
 * @return    The angle in radians in the range [-π, π].
 */
float angle(const Point &p1, const Point &p2);

/**
 * @brief Computes the angle formed by three points with p0 as the origin.
 *
 * @param  p0 The reference point (origin).
 * @param  p1 The first point defining the angle.
 * @param  p2 The second point defining the angle.
 * @return    The angle in radians in the range [-π, π].
 */
float angle(const Point &p0, const Point &p1, const Point &p2);

/**
 * @brief Classify point convexity relative to neighboring curve vertices.
 */
float classify_point(const Point &p_prev,
                     const Point &p,
                     const Point &p_next,
                     const Point &pq);

/**
 * @brief Computes the 2D cross product of vectors (p1 - p0) and (p2 - p0).
 *
 * @param  p0 Common origin point.
 * @param  p1 First destination point.
 * @param  p2 Second destination point.
 * @return    Scalar cross product value.
 */
float cross_product(const Point &p0, const Point &p1, const Point &p2);

/**
 * @brief Computes the 2D cross product of vectors p1 and p2.
 *
 * @param  p1 First vector.
 * @param  p2 Second vector.
 * @return    Scalar cross product (p1.x * p2.y - p1.y * p2.x).
 */
float cross_product(const Point &p1, const Point &p2);

/**
 * @brief Calculates the curvature formed by three points in 2D space.
 *
 * @param  p1 The first point.
 * @param  p2 The second point.
 * @param  p3 The third point.
 * @return    The curvature. Returns 0 if the points are collinear.
 */
float curvature(const Point &p1, const Point &p2, const Point &p3);

/**
 * @brief Calculates the signed curvature formed by three points in 2D space.
 */
float curvature_signed(const Point &p1, const Point &p2, const Point &p3);

/**
 * @brief Calculates the Euclidean distance between two points.
 *
 * @param  p1 The first point.
 * @param  p2 The second point.
 * @return    The distance between the two points.
 */
float distance(const Point &p1, const Point &p2);

/**
 * @brief Calculates the squared Euclidean distance between two points.
 *
 * @param  p1 The first point.
 * @param  p2 The second point.
 * @return    The squared distance between the two points.
 */
float distance_squared(const Point &p1, const Point &p2);

/**
 * @brief Performs a cubic Bezier interpolation.
 *
 * @param  p_start      The first control point (start point).
 * @param  p_ctrl_start The second control point.
 * @param  p_ctrl_end   The third control point.
 * @param  p_end        The fourth control point (end point).
 * @param  t            The interpolation parameter, ranging from 0 to 1.
 * @return              The interpolated point on the Bezier curve.
 */
Point interp_bezier(const Point &p_start,
                    const Point &p_ctrl_start,
                    const Point &p_ctrl_end,
                    const Point &p_end,
                    float        t);

/**
 * @brief Performs a cubic B-spline interpolation.
 *
 * @param  p0 The first control point.
 * @param  p1 The second control point (start of segment).
 * @param  p2 The third control point (end of segment).
 * @param  p3 The fourth control point.
 * @param  t  The interpolation parameter, ranging from 0 to 1.
 * @return    The interpolated point on the B-spline curve.
 */
Point interp_bspline(const Point &p0,
                     const Point &p1,
                     const Point &p2,
                     const Point &p3,
                     float        t);

/**
 * @brief Performs a Catmull-Rom spline interpolation.
 *
 * @param  p0 The first control point.
 * @param  p1 The second control point (start of segment).
 * @param  p2 The third control point (end of segment).
 * @param  p3 The fourth control point.
 * @param  t  The interpolation parameter, ranging from 0 to 1.
 * @return    The interpolated point on the Catmull-Rom spline.
 */
Point interp_catmullrom(const Point &p0,
                        const Point &p1,
                        const Point &p2,
                        const Point &p3,
                        float        t);

/**
 * @brief Performs a De Casteljau algorithm-based interpolation for Bezier
 * curves.
 *
 * @param  points A vector of control points defining the Bezier curve.
 * @param  t      The interpolation parameter, ranging from 0 to 1.
 * @return        The interpolated point on the Bezier curve.
 */
Point interp_decasteljau(const std::vector<Point> &points, float t);

/**
 * @brief Determines the intersection of two bounding boxes.
 *
 * @param  bbox1 The first bounding box defined as `glm::vec4`.
 * @param  bbox2 The second bounding box defined as `glm::vec4`.
 * @return       The intersecting bounding box, or {1, -1, 1, -1} if disjoint.
 */
glm::vec4 intersect_bounding_boxes(const glm::vec4 &bbox1,
                                   const glm::vec4 &bbox2);

/**
 * @brief Checks if a point is within a specified bounding box.
 *
 * @param  p    The point to check.
 * @param  bbox Bounding box {xmin, xmax, ymin, ymax}.
 * @return      true if the point is within the bounding box, false otherwise.
 */
bool is_point_within_bounding_box(const Point &p, const glm::vec4 &bbox);

/**
 * @brief Checks if a point with coordinates (x, y) is within a specified
 * bounding box.
 *
 * @param  x    The x-coordinate.
 * @param  y    The y-coordinate.
 * @param  bbox Bounding box {xmin, xmax, ymin, ymax}.
 * @return      true if the point is within the bounding box, false otherwise.
 */
bool is_point_within_bounding_box(float x, float y, const glm::vec4 &bbox);

/**
 * @brief Linearly interpolates between two points.
 *
 * @param  p1 The starting point.
 * @param  p2 The ending point.
 * @param  t  The interpolation factor (0 <= t <= 1).
 * @return    The interpolated point between `p1` and `p2`.
 */
Point lerp(const Point &p1, const Point &p2, float t);

/**
 * @brief Computes the midpoint displacement with a perpendicular offset.
 *
 * @param  p1             The first point.
 * @param  p2             The second point.
 * @param  orientation    Orientation of the perpendicular displacement.
 * @param  distance_ratio Ratio of displacement relative to segment length.
 * @param  t              Interpolation factor along segment (default 0.5).
 * @return                The displaced midpoint.
 */
Point midpoint(const Point &p1,
               const Point &p2,
               int          orientation,
               float        distance_ratio,
               float        t = 0.5f);

/**
 * @brief Scales a point relative to a center point.
 *
 * @param  p      The input point.
 * @param  scale  Scale factor(s) for x and y axes.
 * @param  center Center of scaling (defaults to (0.5, 0.5)).
 * @return        The scaled point.
 */
Point scale(const Point &p, glm::vec2 scale, glm::vec2 center = {0.5f, 0.5f});

/**
 * @brief Scales a point uniformly relative to a center point.
 */
Point scale(const Point &p, float scale, glm::vec2 center = {0.5f, 0.5f});

/**
 * @brief Computes the intersection point of two 2D segments, if it exists.
 *
 * @param  p1 Start point of first segment.
 * @param  p2 End point of first segment.
 * @param  q1 Start point of second segment.
 * @param  q2 End point of second segment.
 * @return    The intersection point if segments intersect; std::nullopt
 * otherwise.
 */
std::optional<Point> segment_intersection(const Point &p1,
                                          const Point &p2,
                                          const Point &q1,
                                          const Point &q2);

/**
 * @brief Determines the relative side of a query point with respect to a curve
 * segment.
 *
 * @param  p1      Preceding curve point.
 * @param  p2      Reference curve point.
 * @param  p3      Succeeding curve point.
 * @param  p_query Query point to test.
 * @return         +1.0f (left/CCW), -1.0f (right/CW), 0.0f (on tangent).
 */
float side(const Point &p1,
           const Point &p2,
           const Point &p3,
           const Point &p_query);

/**
 * @brief Sorts a vector of points in ascending order (x, then y, then v).
 *
 * @param points Vector of points to sort in-place.
 */
void sort_points(std::vector<Point> &points);

/**
 * @brief Calculates the unsigned area of a triangle formed by three points.
 *
 * @param  p1 The first point.
 * @param  p2 The second point.
 * @param  p3 The third point.
 * @return    The absolute area.
 */
float triangle_area(const Point &p1, const Point &p2, const Point &p3);

/**
 * @brief Calculates the signed area of a triangle formed by three points.
 *
 * @param  p1 The first point.
 * @param  p2 The second point.
 * @param  p3 The third point.
 * @return    The signed area (positive = CCW, negative = CW).
 */
float triangle_area_signed(const Point &p1, const Point &p2, const Point &p3);

/**
 * @brief Constructs a 4D bounding box for a unit square {0, 1, 0, 1}.
 *
 * @return glm::vec4 Unit square bounding box.
 */
glm::vec4 unit_square_bbox();

} // namespace hmap