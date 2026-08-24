/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file array.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Absolute value.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array abs(const Array &array);

/**
 * @brief Smooth absolute value.
 *
 * @param  array  Input array.
 * @param  mu     Smoothing coefficient.
 * @param  vshift Reference value for the zero of the absolute value.
 * @return        Output array.
 *
 * **Example**
 * @include ex_abs_smooth.cpp
 *
 * **Result**
 * @image html ex_abs_smooth.png
 */
Array abs_smooth(const Array &array, float mu, const Array &vshift);
Array abs_smooth(const Array &array, float mu, float vshift); ///< @overload
Array abs_smooth(const Array &array, float mu);               ///< @overload

/**
 * @brief Almost unit identity function.
 *
 * Maps [0, 1] to itself with zero derivative at 0 and unit derivative at 1. See
 * Inigo Quilez: https://iquilezles.org/articles/functions/
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array almost_unit_identity(const Array &array);

/**
 * @brief Almost unit identity function (C2 variant).
 *
 * Ensures smoother behavior by also enforcing second-order continuity at x = 1.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array almost_unit_identity_c2(const Array &array);

/**
 * @brief Arctangent.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array atan(const Array &array);

/**
 * @brief Element-wise arctangent with quadrant handling.
 *
 * Computes atan2(y, x) element-wise.
 *
 * @param  y Numerator array.
 * @param  x Denominator array.
 * @return   Output array.
 */
Array atan2(const Array &y, const Array &array);

/**
 * @brief Cosine.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array cos(const Array &array);

/**
 * @brief Exponential.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array exp(const Array &array);

/**
 * @brief Gaussian decay.
 *
 * @param  array Input array.
 * @param  sigma Gaussian half-width.
 * @return       Output array.
 *
 * **Example**
 * @include ex_gaussian_decay.cpp
 *
 * **Result**
 * @image html ex_gaussian_decay.png
 */
Array gaussian_decay(const Array &array, float sigma);

/**
 * @brief Hypotenuse.
 *
 * Computes sqrt(array1² + array2²) element-wise.
 *
 * @param  array1 First array.
 * @param  array2 Second array.
 * @return        Output array.
 */
Array hypot(const Array &array1, const Array &array2);

/**
 * @brief Binary mask of elements equal to a value.
 *
 * @param  array Input array.
 * @param  value Reference value.
 * @return       Binary mask array.
 */
Array is_equal(const Array &array, float value);

/**
 * @brief Binary mask of non-zero elements.
 *
 * @param  array Input array.
 * @return       Binary mask array.
 */
Array is_non_zero(const Array &array);

/**
 * @brief Binary mask of zero elements.
 *
 * @param  array Input array.
 * @return       Binary mask array.
 */
Array is_zero(const Array &array);

/**
 * @brief Linear interpolation.
 *
 * @param  array1 First array.
 * @param  array2 Second array.
 * @param  t      Interpolation factor in [0,1].
 * @return        Interpolated array.
 */
Array lerp(const Array &array1, const Array &array2, const Array &t);
Array lerp(const Array &array1, const Array &array2, float t); ///< @overload

/**
 * @brief Base-10 logarithm.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array log10(const Array &array);

/**
 * @brief Power function.
 *
 * @param  array Input array.
 * @param  exp   Exponent.
 * @return       Output array.
 */
Array pow(const Array &array, float exp);

/**
 * @brief Generalized Rvachev R-function for smooth minimum / intersection
 * (exact zero-set).
 *
 * @param  array1 First array.
 * @param  array2 Second array.
 * @param  alpha  Parameter in (-1, 1] controlling blending sharpness (default
 *                0.0).
 * @return        Output array.
 */
Array r_min(const Array &array1, const Array &array2, float alpha = 0.f);

/**
 * @brief Generalized Rvachev R-function for smooth maximum / union (exact
 * zero-set).
 *
 * @param  array1 First array.
 * @param  array2 Second array.
 * @param  alpha  Parameter in (-1, 1] controlling blending sharpness (default
 *                0.0).
 * @return        Output array.
 */
Array r_max(const Array &array1, const Array &array2, float alpha = 0.f);

/**
 * @brief Sigmoid function.
 *
 * Applies the generalized sigmoid transformation element-wise:
 * \f[
 * y = v_{\min} + \frac{v_{\max} - v_{\min}}
 *     {1 + \exp\left(-\frac{x - x_0}{\text{width}}\right)}
 * \f]
 *
 * @param  array Input array.
 * @param  width Controls steepness.
 * @param  vmin  Minimum output value.
 * @param  vmax  Maximum output value.
 * @param  x0    Center offset.
 * @return       Output array.
 *
 * **Example**
 * @include ex_sigmoid.cpp
 *
 * **Result**
 * @image html ex_sigmoid.png
 */
Array sigmoid(const Array &array,
              float        width = 1.f,
              float        vmin = 0.f,
              float        vmax = 1.f,
              float        x0 = 0.f);

/**
 * @brief Sine.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array sin(const Array &array);

/**
 * @brief Cubic smoothstep.
 *
 * @param  array Input array.
 * @param  vmin  Lower bound.
 * @param  vmax  Upper bound.
 * @return       Output array.
 *
 * **Example**
 * @include ex_smoothstep.cpp
 *
 * **Result**
 * @image html ex_smoothstep.png
 */
Array smoothstep3(const Array &array, float vmin = 0.f, float vmax = 1.f);

/**
 * @brief Cubic smoothstep (zero derivative at 0 only).
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array smoothstep3_lower(const Array &array); ///< @overload

/**
 * @brief Cubic smoothstep (zero derivative at 1 only).
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array smoothstep3_upper(const Array &array); ///< @overload

/**
 * @brief Quintic smoothstep.
 *
 * @param  array Input array.
 * @param  vmin  Lower bound.
 * @param  vmax  Upper bound.
 * @return       Output array.
 *
 * **Example**
 * @include ex_smoothstep.cpp
 *
 * **Result**
 * @image html ex_smoothstep.png
 */
Array smoothstep5(const Array &array, float vmin = 0.f, float vmax = 1.f);
Array smoothstep5(const Array &array, const Array &vmin, const Array &vmax);

/**
 * @brief Quintic smoothstep (zero derivative at 0 only).
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array smoothstep5_lower(const Array &array);

/**
 * @brief Quintic smoothstep (zero derivative at 1 only).
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array smoothstep5_upper(const Array &array);

/**
 * @brief Seventh-order smoothstep.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array smoothstep7(const Array &array);

/**
 * @brief Square root.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array sqrt(const Array &array);

/**
 * @brief Safe square root.
 *
 * Returns sqrt(v) for positive values and 0 otherwise.
 *
 * @param  array Input array.
 * @return       Output array.
 */
Array sqrt_safe(const Array &array);

/**
 * @brief Linear threshold.
 *
 * @param  array Input array.
 * @param  x0    Lower threshold.
 * @param  x1    Upper threshold.
 * @return       Output array.
 */
Array threshold(const Array &array, float x0, float x1);

/**
 * @brief Smooth threshold.
 *
 * @param  array Input array.
 * @param  x0    Lower threshold.
 * @param  x1    Upper threshold.
 * @return       Output array.
 */
Array threshold_smooth(const Array &array, float x0, float x1);

} // namespace hmap