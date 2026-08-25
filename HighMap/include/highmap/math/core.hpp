/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file core.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2023
 */
#pragma once
#include <cmath>
#include <cstdint>

namespace hmap
{

/**
 * @brief Smooth absolute value.
 *
 * @param a  Input value.
 * @param mu Smoothing coefficient.
 */
float abs_smooth(float a, float mu);

/**
 * @brief Almost unit identity function.
 *
 * Maps [0, 1] to itself with zero derivative at 0 and unit derivative at 1. See
 * Inigo Quilez: https://iquilezles.org/articles/functions/
 *
 * @param x Input value.
 */
float almost_unit_identity(float x);

/**
 * @brief Almost unit identity function (C2 variant).
 *
 * Ensures smoother behavior by also enforcing second-order continuity at x = 1.
 *
 * @param x Input value.
 */
float almost_unit_identity_c2(float x);

/**
 * @brief Approximate hypotenuse.
 *
 * @param  a First value.
 * @param  b Second value.
 * @return   Approximate sqrt(a² + b²).
 */
inline float approx_hypot(float a, float b)
{
  a = std::abs(a);
  b = std::abs(b);
  if (a > b) std::swap(a, b);
  return 0.414f * a + b;
}

/**
 * @brief Fast inverse square root approximation.
 *
 * @param  a Input value.
 * @return   Approximate 1 / sqrt(a).
 */
inline float approx_rsqrt(float a)
{
  union
  {
    float    f;
    uint32_t i;
  } conv = {.f = a};

  conv.i = 0x5f3759df - (conv.i >> 1);
  conv.f *= 1.5F - (a * 0.5F * conv.f * conv.f);

  return conv.f;
}

/**
 * @brief Integer ceiling division.
 *
 * @param  a Numerator (>= 0).
 * @param  b Denominator (> 0).
 * @return   ceil(a / b).
 */
int ceil_div(int a, int b);

/**
 * @brief Fast exponential approximation.
 *
 * @param  x Input value.
 * @return   Approximate exp(x).
 */
float fast_exp(float x);

/**
 * @brief Fast natural logarithm approximation.
 *
 * @param  x Input value (must be > 0).
 * @return   Approximate ln(x).
 */
float fast_log(float x);

/**
 * @brief Gain curve remapping in [0,1].
 *
 * @param x      Input value in [0,1].
 * @param factor Contrast/gain factor.
 */
float gain(float x, float factor);

/**
 * @brief Highest power of 2 less than or equal to n.
 *
 * @param  n Input value.
 * @return   Power of 2.
 */
int highest_power_of_2(int n);

/**
 * @brief Linear interpolation.
 *
 * @param a First value.
 * @param b Second value.
 * @param t Interpolation factor in [0,1].
 */
float lerp(float a, float b, float t);

/**
 * @brief Smooth asymmetric bell-shaped curve on [0,1].
 *
 * Parameters @p a and @p b control left and right curvature. The function is
 * @param x Input value.
 * @param a Left curvature parameter.
 * @param b Right curvature parameter.
 */
float power_curve(float x, float a, float b);

/**
 * @brief Generalized Rvachev R-function for smooth minimum / intersection
 * (exact zero-set).
 *
 * \f[
 * R_{\min}(f_1, f_2, \alpha) = \frac{f_1 + f_2 - \sqrt{f_1^2 + f_2^2 - 2\alpha
 * f_1 f_2}}{1 + \alpha} \f]
 *
 * @param  f1    First value.
 * @param  f2    Second value.
 * @param  alpha Parameter in (-1, 1] controlling blending sharpness (default
 *               0.0). alpha = 0 gives standard Euclidean R-min, alpha = 1 gives
 *               exact min.
 * @return       Smooth minimum value.
 */
float r_min(float f1, float f2, float alpha = 0.f);

/**
 * @brief Generalized Rvachev R-function for smooth maximum / union (exact
 * zero-set).
 *
 * \f[
 * R_{\max}(f_1, f_2, \alpha) = \frac{f_1 + f_2 + \sqrt{f_1^2 + f_2^2 - 2\alpha
 * f_1 f_2}}{1 + \alpha} \f]
 *
 * @param  f1    First value.
 * @param  f2    Second value.
 * @param  alpha Parameter in (-1, 1] controlling blending sharpness (default
 *               0.0). alpha = 0 gives standard Euclidean R-max, alpha = 1 gives
 *               exact max.
 * @return       Smooth maximum value.
 */
float r_max(float f1, float f2, float alpha = 0.f);

/**
 * @brief Sigmoid function.
 *
 * Generalized logistic function:
 * \f[
 * y = v_{\min} + \frac{v_{\max} - v_{\min}}
 *     {1 + \exp\left(-\frac{x - x_0}{\text{width}}\right)}
 * \f]
 *
 * @param  x     Input value.
 * @param  width Controls steepness.
 * @param  vmin  Minimum output value.
 * @param  vmax  Maximum output value.
 * @param  x0    Center offset.
 * @return       Sigmoid result in [vmin, vmax].
 */
float sigmoid(float x,
              float width = 1.f,
              float vmin = 0.f,
              float vmax = 1.f,
              float x0 = 0.f);

/**
 * @brief Cubic smoothstep function.
 *
 * @param x Input value.
 */
float smoothstep3(float x);

/**
 * @brief Cubic smoothstep (zero derivative at 0 only).
 *
 * @param x Input value.
 */
float smoothstep3_lower(float x);

/**
 * @brief Cubic smoothstep (zero derivative at 1 only).
 *
 * @param x Input value.
 */
float smoothstep3_upper(float x);

/**
 * @brief Quintic smoothstep function.
 *
 * @param x Input value.
 */
float smoothstep5(float x);

/**
 * @brief Quintic smoothstep (zero derivative at 0 only).
 *
 * @param x Input value.
 */
float smoothstep5_lower(float x);

/**
 * @brief Quintic smoothstep (zero derivative at 1 only).
 *
 * @param x Input value.
 */
float smoothstep5_upper(float x);

/**
 * @brief Seventh-order smoothstep function.
 *
 * @param x Input value.
 */
float smoothstep7(float x);

/**
 * @brief Linear threshold.
 *
 * @param  x  Input value.
 * @param  x0 Lower bound.
 * @param  x1 Upper bound.
 * @return    0 if x < x0, 1 if x > x1, linear otherwise.
 */
float threshold(float x, float x0, float x1);

/**
 * @brief Smooth threshold using smoothstep.
 *
 * @param x  Input value.
 * @param x0 Lower bound.
 * @param x1 Upper bound.
 */
float threshold_smooth(float x, float x0, float x1);

/**
 * @brief Computes a trapezoidal pulse function.
 *
 * Returns a value in [0, 1] with linear ramps of size `width`
 * around the interval [x0, x1].
 *
 * @param x     Input value.
 * @param x0    Start of the plateau.
 * @param x1    End of the plateau.
 * @param width Width of the rising and falling ramps.
 */
float trapeze(float x, float x0, float x1, float width);

/**
 * @brief Smooth trapezoidal pulse using cubic smoothstep interpolation.
 *
 * @param x     Input value.
 * @param x0    Start of the plateau.
 * @param x1    End of the plateau.
 * @param width Width of the rising and falling ramps.
 */
float trapeze_smooth(float x, float x0, float x1, float width);

/**
 * @brief Triangle function between bounds.
 *
 * @param  x    Input value.
 * @param  vmin Lower bound.
 * @param  vmax Upper bound.
 * @return      0 outside bounds, peak 1 at midpoint.
 */
float triangle(float x, float vmin, float vmax);

} // namespace hmap