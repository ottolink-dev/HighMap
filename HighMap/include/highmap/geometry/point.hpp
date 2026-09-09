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

  /**
   * @brief Default constructor initializing the point to (0, 0, 0).
   */
  Point() : x(0.f), y(0.f), v(0.f)
  {
  }

  /**
   * @brief Parameterized constructor initializing the point to given values.
   * @param x The x-coordinate of the point.
   * @param y The y-coordinate of the point.
   * @param v The value at the point.
   */
  Point(float x, float y, float v = 0.f) : x(x), y(y), v(v)
  {
  }

  /**
   * @brief Equality operator to check if two points are the same.
   *
   * @param  other The point to compare with.
   * @return       true if the points are equal, false otherwise.
   */
  bool operator==(const Point &other) const
  {
    return (x == other.x && y == other.y && v == other.v);
  }

  /**
   * @brief Inequality operator to check if two points are different.
   *
   * @param  other The point to compare with.
   * @return       true if the points are not equal, false otherwise.
   */
  bool operator!=(const Point &other) const
  {
    return !(*this == other);
  }

  /**
   * @brief Adds two points.
   * @param  other The point to add.
   * @return       The result of adding the two points.
   */
  Point operator+(const Point &other) const
  {
    return Point(x + other.x, y + other.y, v + other.v);
  }

  /**
   * @brief Subtracts two points.
   * @param  other The point to subtract.
   * @return       The result of subtracting the other point from this point.
   */
  Point operator-(const Point &other) const
  {
    return Point(x - other.x, y - other.y, v - other.v);
  }

  /**
   * @brief Multiplies the point by a scalar.
   * @param  scalar The scalar to multiply by.
   * @return        The result of the multiplication.
   */
  Point operator*(float scalar) const
  {
    return Point(x * scalar, y * scalar, v * scalar);
  }

  /**
   * @brief Divides the point by a scalar.
   * @param  scalar The scalar to divide by.
   * @return        The result of the division.
   */
  Point operator/(float scalar) const
  {
    return Point(x / scalar, y / scalar, v / scalar);
  }

  /**
   * @brief Scalar multiplication (scalar * Vec2).
   *
   * Multiplies each component of the vector by a scalar value. This function
   * allows expressions where the scalar is on the left side of the
   * multiplication operator.
   *
   * @param  scalar The scalar value to multiply with.
   * @param  point  The vector to multiply.
   * @return        Point A new vector with each component multiplied by the
   *                scalar.
   */
  friend Point operator*(float scalar, const Point &point)
  {
    return Point(scalar * point.x, scalar * point.y, scalar * point.v);
  }

  /**
   * @brief Prints the coordinates and value of the Point object.
   *
   * This function outputs the Point's x, y coordinates, and an additional value
   * `v` to the standard output in the format `(x, y, v)`, followed by a
   * newline.
   */
  void print();

  /**
   * @brief Updates the point's value based on bilinear interpolation from an
   * array.
   *
   * This function updates the value of the `Point` object by performing
   * bilinear interpolation on the input `Array` using the coordinates of the
   * `Point` and a given bounding box. The point's coordinates are first
   * normalized to the unit interval using the provided bounding box. Then,
   * these normalized coordinates are scaled to the array's dimensions and used
   * to fetch the value from the array through bilinear interpolation.
   *
   * If the normalized coordinates fall outside the bounds of the array, the
   * point's value is set to zero.
   *
   * @param array The input `Array` from which the value is interpolated. The
   * `Array` should support bilinear interpolation.
   * @param bbox  Bounding box used for normalizing the `Point`'s coordinates.
   *              This box is defined by a `glm::vec4` containing minimum and
   *              maximum values for both x and y dimensions in the format
   * `{xmin, xmax, ymin, ymax}`.
   *
   * @note If the coordinates are outside the bounds of the array after scaling,
   * the point's value is set to zero.
   */
  void set_value_from_array(const Array &array, glm::vec4 bbox);
};

/**
 * @brief Computes the angle between two points relative to the x-axis.
 *
 * This function calculates the angle formed by the vector from `p1` to `p2`
 * with respect to the x-axis. The angle is measured in radians and is oriented
 * in the counter-clockwise direction.
 *
 * @param  p1 The starting point of the vector.
 * @param  p2 The ending point of the vector.
 * @return    The angle in radians between the vector formed by `p1` to `p2` and
 *            the x-axis. The angle is in the range [-π, π].
 */
float angle(const Point &p1, const Point &p2);

/**
 * @brief Computes the angle formed by three points with the reference point as
 * the origin.
 *
 * Given three points \( p0 \), \( p1 \), and \( p2 \), this function calculates
 * the angle formed between the vectors \( p0 \rightarrow p2 \) and \( p0
 * \rightarrow p1 \). The angle is oriented in the 2D plane and measured in
 * radians.
 *
 * @param  p0 The reference point (origin of the angle measurement).
 * @param  p1 The first point defining the angle.
 * @param  p2 The second point defining the angle.
 * @return    The angle between the vectors \( p0 \rightarrow p2 \) and \( p0
 * \rightarrow p1 \) in radians.
 *
 * @note The angle is calculated using the 2D vectors defined by \( p0
 * \rightarrow p1 \) and \( p0 \rightarrow p2 \). The result will be in the
 * range \([-π, π]\).
 */
float angle(const Point &p0, const Point &p1, const Point &p2);

float classify_point(const Point &p_prev,
                     const Point &p,
                     const Point &p_next,
                     const Point &pq);

/**
 * @brief Computes the 2D cross product of vectors formed by three points.
 *
 * Given three points p0, p1, and p2 in 2D space, this function computes the
 * scalar cross product of the vectors (p1 - p0) and (p2 - p0). In 2D, the cross
 * product is a scalar value and is often used to determine the orientation of
 * the three points or the signed area of the parallelogram they form.
 *
 * @param  p0 The first point, which serves as the common origin of the two
 *            vectors.
 * @param  p1 The second point, forming the first vector p1 - p0.
 * @param  p2 The third point, forming the second vector p2 - p0.
 * @return    The scalar value of the 2D cross product.
 *         - Positive if the points p0, p1, p2 are oriented counterclockwise.
 *         - Negative if they are oriented clockwise.
 *         - Zero if the points are collinear.
 */
float cross_product(const Point &p0, const Point &p1, const Point &p2);
float cross_product(const Point &p1, const Point &p2);

/**
 * @brief Calculates the curvature formed by three points in 2D space.
 *
 * @param  p1 The first point of the triangle.
 * @param  p2 The second point of the triangle.
 * @param  p3 The third point of the triangle.
 * @return    The curvature. Returns 0 if the points are collinear.
 */
float curvature(const Point &p1, const Point &p2, const Point &p3);

/* see hmap::curvature */
float curvature_signed(const Point &p1, const Point &p2, const Point &p3);

/**
 * @brief Calculates the distance between two points.
 *
 * This function computes the Euclidean distance between two points in 2D space.
 *
 * @param  p1 The first point.
 * @param  p2 The second point.
 * @return    The distance between the two points.
 */
float distance(const Point &p1, const Point &p2);

/**
 * @brief Performs a cubic Bezier interpolation.
 *
 * Interpolates a point on a cubic Bezier curve defined by the control points
 * `p_start`, `p_ctrl_start`, `p_ctrl_end`, and `p_end`, using the parameter
 * `t`. The parameter `t` should be within the range [0, 1], where `t = 0`
 * corresponds to the start point `p_start` and `t = 1` corresponds to the end
 * point `p_end`.
 *
 * @param  p_start      The first control point (start point).
 * @param  p_ctrl_start The second control point.
 * @param  p_ctrl_end   The third control point.
 * @param  p_end        The fourth control point (end point).
 * @param  t            The interpolation parameter, ranging from 0 to 1.
 * @return              The interpolated point on the Bezier curve.
 *
 * **Example**
 * @include ex_point_interp.cpp
 *
 * **Result**
 * @image html ex_point_interp.png
 */
Point interp_bezier(const Point &p_start,
                    const Point &p_ctrl_start,
                    const Point &p_ctrl_end,
                    const Point &p_end,
                    float        t);

/**
 * @brief Performs a cubic B-spline interpolation.
 *
 * Interpolates a point on a cubic B-spline curve using the control points
 * `p0`, `p1`, `p2`, and `p3`. The points `p1` and `p2` define the segment of
 * the curve to be interpolated, while `p0` and `p3` are used as additional
 * control points. The parameter `t` should be within the range [0, 1], where
 * `t = 0` corresponds to the start point `p1` and `t = 1` corresponds to the
 * end point `p2`.
 *
 * @param  p0 The first control point (influence for the start point).
 * @param  p1 The second control point (start point of the segment).
 * @param  p2 The third control point (end point of the segment).
 * @param  p3 The fourth control point (influence for the end point).
 * @param  t  The interpolation parameter, ranging from 0 to 1.
 * @return    The interpolated point on the B-spline curve.
 *
 * **Example**
 * @include ex_point_interp.cpp
 *
 * **Result**
 * @image html ex_point_interp.png
 */
Point interp_bspline(const Point &p0,
                     const Point &p1,
                     const Point &p2,
                     const Point &p3,
                     float        t);

/**
 * @brief Performs a Catmull-Rom spline interpolation.
 *
 * Interpolates a point on a Catmull-Rom spline defined by the points `p0`,
 * `p1`, `p2`, and `p3`. The points `p1` and `p2` define the segment of the
 * curve to be interpolated, while `p0` and `p3` are used as additional control
 * points. The parameter `t` should be within the range [0, 1], where `t = 0`
 * corresponds to the start point `p1` and `t = 1` corresponds to the end point
 * `p2`.
 *
 * @param  p0 The first control point (influence for the start point).
 * @param  p1 The second control point (start point of the segment).
 * @param  p2 The third control point (end point of the segment).
 * @param  p3 The fourth control point (influence for the end point).
 * @param  t  The interpolation parameter, ranging from 0 to 1.
 * @return    The interpolated point on the Catmull-Rom spline.
 *
 * **Example**
 * @include ex_point_interp.cpp
 *
 * **Result**
 * @image html ex_point_interp.png
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
 * Interpolates a point on a Bézier curve defined by a set of control points
 * using De Casteljau's algorithm. The algorithm is recursive and provides a
 * stable and numerically robust method to evaluate Bézier curves at a given
 * parameter `t`.
 *
 * The parameter `t` should be within the range [0, 1], where `t = 0`
 * corresponds to the first control point and `t = 1` corresponds to the last
 * control point.
 *
 * @param  points A vector of control points defining the Bézier curve.
 * @param  t      The interpolation parameter, ranging from 0 to 1.
 * @return        The interpolated point on the Bézier curve.
 *
 * **Example**
 * @include ex_point_interp.cpp
 *
 * **Result**
 * @image html ex_point_interp.png
 */
Point interp_decasteljau(const std::vector<Point> &points, float t);

/**
 * @brief Determines the intersection of two bounding boxes.
 *
 * This function calculates the intersection of two bounding boxes, `bbox1` and
 * `bbox2`. If they intersect, it returns the intersecting bounding box. If they
 * are disjoint, it returns `std::nullopt`.
 *
 * @param  bbox1 The first bounding box defined as `glm::vec4`.
 * @param  bbox2 The second bounding box defined as `glm::vec4`.
 * @return       An `std::optional<glm::vec4>` containing the intersecting
 *               bounding box if an intersection exists; `std::nullopt`
 *               otherwise.
 */
glm::vec4 intersect_bounding_boxes(const glm::vec4 &bbox1,
                                   const glm::vec4 &bbox2);

/**
 * @brief Checks if a point is within a specified bounding box.
 *
 * This function determines if a given point `p` lies within the rectangular
 * bounding box defined by `bbox`.
 *
 * @param  p    The point to check, represented as a `Point` with `x` and `y`
 * coordinates.
 * @param  bbox The bounding box defined as a `glm::vec4`, where `a` and `b`
 * are the horizontal boundaries (min and max x), and `c` and `d` are the
 * vertical boundaries (min and max y).
 * @return      `true` if the point is within the bounding box; `false`
 *              otherwise.
 */
bool is_point_within_bounding_box(Point p, glm::vec4 bbox);

/**
 * @brief Checks if a point is within a specified bounding box.
 *
 * This function determines if a point with coordinates `(x, y)` lies within the
 * rectangular bounding box defined by `bbox`.
 *
 * @param  x    The x-coordinate of the point to check.
 * @param  y    The y-coordinate of the point to check.
 * @param  bbox The bounding box defined as a `glm::vec4`, where `a` and `b`
 * are the horizontal boundaries (min and max x), and `c` and `d` are the
 * vertical boundaries (min and max y).
 * @return      `true` if the point is within the bounding box; `false`
 *              otherwise.
 */
bool is_point_within_bounding_box(float x, float y, glm::vec4 bbox);

/**
 * @brief Linearly interpolates between two points.
 *
 * This function performs linear interpolation between two points based on a
 * given factor. The interpolation factor `t` should be in the range [0, 1],
 * where 0 corresponds to the first point and 1 corresponds to the second point.
 *
 * @param  p1 The starting point.
 * @param  p2 The ending point.
 * @param  t  The interpolation factor (0 <= t <= 1).
 * @return    The interpolated point between `p1` and `p2`.
 *
 * @note If `t` is outside the range [0, 1], the function will still return a
 * point outside the segment defined by `p1` and `p2`.
 */
Point lerp(const Point &p1, const Point &p2, float t);

/**
 * @brief Computes the midpoint displacement in 1D with a perpendicular
 * displacement.
 *
 * This function generates a midpoint between two points `p1` and `p2` based on
 * a linear interpolation parameter `t`. It then displaces this midpoint in the
 * direction perpendicular to the line segment formed by `p1` and `p2`. The
 * displacement distance is determined by the `distance_ratio` parameter and the
 * `orientation` (which can be -1 or 1).
 *
 * @param  p1             The first point.
 * @param  p2             The second point.
 * @param  orientation    Determines the direction of the perpendicular
 *                        displacement.
 *                    - A value of `1` displaces the midpoint in the positive
 * perpendicular direction.
 *                    - A value of `-1` displaces the midpoint in the negative
 * perpendicular direction.
 * @param  distance_ratio The ratio of the displacement distance relative to the
 *                        length of the line segment `p1p2`.
 * @param  t              The interpolation factor (default is 0.5) for
 *                        computing the midpoint between `p1` and `p2`, before
 *                        the displacement.
 * @return                A `Point` representing the displaced midpoint.
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
 * @return        Point  The scaled point.
 */
Point scale(const Point &p, glm::vec2 scale, glm::vec2 center = {0.5f, 0.5f});

Point scale(const Point &p, float scale, glm::vec2 center = {0.5f, 0.5f});

/**
 * @brief Computes the intersection point of two 2D segments, if it exists.
 *
 * @param  p1 Start point of the first segment.
 * @param  p2 End point of the first segment.
 * @param  q1 Start point of the second segment.
 * @param  q2 End point of the second segment.
 * @return    std::optional<Point> The intersection point if the segments
 *            intersect; std::nullopt otherwise.
 */
std::optional<Point> segment_intersection(const Point &p1,
                                          const Point &p2,
                                          const Point &q1,
                                          const Point &q2);

/**
 * @brief Determines the relative side of a query point with respect to a curve
 * segment at a given point.
 *
 * @param  p1 A point preceding @p p2 on the curve.
 * @param  p2 The reference point on the curve where the tangent is evaluated.
 * @param  p3 A point following @p p2 on the curve.
 * @param  pq The query point to test.
 *
 * @return    float
 * - +1.0f if the query point lies to the left of the curve (counter-clockwise
 * side),
 * - -1.0f if the query point lies to the right of the curve (clockwise side),
 * -  0.0f if the query point lies exactly on the tangent line at @p p2.
 */
float side(const Point &p1,
           const Point &p2,
           const Point &p3,
           const Point &p_query);

/**
 * @brief Sorts a vector of points in ascending order based on their
 * coordinates.
 *
 * This function sorts the provided vector of `Point` objects first by the
 * x-coordinate and then by the y-coordinate in ascending order. The sorting is
 * done in-place.
 *
 * @param points A vector of `Point` objects to be sorted. The vector is
 *               modified directly with the points arranged in ascending order.
 */
void sort_points(std::vector<Point> &points);

/**
 * @brief Calculates the area of a triangle formed by three points in 2D space.
 *
 * This function uses the determinant method to compute the absolute area of a
 * triangle given three points (p1, p2, p3). It assumes the points are specified
 * as 2D coordinates.
 *
 * @param  p1 The first point of the triangle.
 * @param  p2 The second point of the triangle.
 * @param  p3 The third point of the triangle.
 * @return    The area of the triangle as a floating-point value. Returns 0 if
 *            the points are collinear.
 */
float triangle_area(const Point &p1, const Point &p2, const Point &p3);

/* see hmap::triangle_area */
float triangle_area_signed(const Point &p1, const Point &p2, const Point &p3);

/**
 * @brief Constructs a 4D bounding box for a unit square.
 *
 * This function returns a glm::vec4 object representing the bounding box of a
 * unit square, with components (min_x, max_x, min_y, max_y) set to (0.f, 1.f,
 * 0.f, 1.f).
 *
 * @return glm::vec4 A 4D vector representing the unit square's bounding box.
 */
glm::vec4 unit_square_bbox();

} // namespace hmap