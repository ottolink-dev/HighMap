/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file interpolator_curve.hpp
 * @author * Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for the `InterpolatorCurve` class, which provides
 * functionality for curve interpolation using various methods such as Bezier,
 * B-spline, and Catmull-Rom.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <functional>

#include "highmap/geometry/point.hpp"

namespace hmap
{
/**
 * @brief Enumeration for specifying the interpolation method for curves.
 *
 * Defines the available methods for curve interpolation:
 * - `BEZIER`: Requires the number of points to be a multiple of 4 minus 1, like
 * (4 * n - 1).
 * - `BSPLINE`: B-spline interpolation method.
 * - `CATMULLROM`: Catmull-Rom spline interpolation method.
 * - `DECASTELJAU`: Uses De Casteljau's algorithm to compute Bézier curves at a
 * given parameter. This is an alternative to the standard Bézier curve method.
 * - `POINTS_LERP`: Linear interpolation between points.
 */
enum InterpolationMethodCurve : int
{
  BEZIER,      ///< Bezier curve interpolation
  BSPLINE,     ///< B-spline interpolation
  CATMULLROM,  ///< Catmull-Rom spline interpolation
  DECASTELJAU, ///< De Casteljau algorithm for Bézier curve computation
  POINTS_LERP  ///< Linear interpolation between points
};

/**
 * @brief Class for performing curve interpolation on a set of points.
 *
 * The `InterpolatorCurve` class allows interpolation over a set of points using
 * a specified interpolation method (Bezier, B-spline, Catmull-Rom, or linear
 * interpolation).
 */
class InterpolatorCurve
{
public:
  /**
   * @brief Constructs an `InterpolatorCurve` object.
   *
   * Initializes the interpolator with a vector of points and an interpolation
   * method.
   *
   * @param points A vector of `Point` objects to interpolate.
   * @param method The interpolation method to use. Default is `POINTS_LERP`.
   */
  InterpolatorCurve(
      std::vector<Point>       points,
      InterpolationMethodCurve method = InterpolationMethodCurve::POINTS_LERP);

  /**
   * @brief Performs interpolation over the provided parameter values.
   *
   * Given a vector of parameter values `t`, this function returns a vector of
   * interpolated points corresponding to those parameter values.
   *
   * @param  t A vector of parameter values, typically in the range [0, 1].
   * @return   A vector of interpolated `Point` objects.
   */
  std::vector<Point> operator()(std::vector<float> t) const;

private:
  std::vector<Point>       points_data; ///< The set of points to interpolate.
  InterpolationMethodCurve method; ///< The interpolation method being used.
  std::vector<float>
      arc_length; ///< Normalized cumulative distance between points.
  std::function<Point(float)>
      interp; ///< Function to perform interpolation based on `method`.

  /**
   * @brief Computes segment interpolation parameters for a given parameter
   * value.
   *
   * This function calculates the appropriate segment index and local parameter
   * value within that segment for a given global parameter `t`.
   *
   * @param t     The global parameter value, typically in the range [0, 1].
   * @param ileft Reference to store the index of the left segment point.
   * @param u     Reference to store the local parameter value within the
   *              segment.
   */
  void get_segment_interpolation_parameters(float t, size_t &ileft, float &u);
};

} // namespace hmap
