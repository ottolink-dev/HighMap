/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file statistics.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2026
 */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{
/**
 * @brief Normalization methods.
 */
enum NormalizationMethod : int
{
  NM_MIN_MAX,     ///<
  NM_STANDARDIZE, ///<
  NM_ROBUST,      ///<
};

/**
 * @brief Estimate the isotropic dominant length scale of a 2D heightmap.
 *
 * Computes the e-folding lag of the autocorrelation function (averaged over i
 * and j directions) via binary search. Returns `max(ni, nj)` for flat arrays.
 *
 * @param  array            Input 2D heightmap, accessed via array(i, j).
 * @param  max_lag_fraction Search range as a fraction of the shorter dimension.
 * @return                  Length scale in grid cells.
 *
 * **Example**
 * @include ex_autocorr_length_scale.cpp
 *
 * See unit tests: @ref test_autocorr_length_scale.cpp
 */
float autocorr_length_scale(const Array &array, float max_lag_fraction = 0.4f);

/**
 * @brief Estimate axis-aligned dominant length scales of a 2D heightmap.
 *
 * Independently computes the e-folding autocorrelation lag along the i and j
 * grid axes. Use when the field is anisotropic but aligned with the grid.
 * Returns `{ni, nj}` for flat arrays.
 *
 * @param  array            Input 2D heightmap, accessed via array(i, j).
 * @param  max_lag_fraction Search range as a fraction of each axis length.
 * @return                  {scale_i, scale_j} length scales in grid cells.
 *
 * See unit tests: @ref test_autocorr_length_scale.cpp
 */
glm::vec2 autocorr_length_scale_axial(const Array &array,
                                      float        max_lag_fraction = 0.4f);

/**
 * @brief Apply linear regression for detrending of a 2D array.
 *
 * This function performs detrending on the input array by applying linear
 * regression separately to each row and column, removing trends from the data.
 *
 * @param  array Input 2D array to be detrended.
 * @return       Array The detrended output array.
 *
 * **Example**
 * @include ex_detrend.cpp
 *
 * **Result**
 * @image html ex_detrend.png
 */
Array detrend_reg(const Array &array);

/**
 * @brief In-place normalization of an array.
 *
 * Applies a normalization transform to the input array depending on the
 * selected method. The operation modifies the array directly.
 *
 * Supported methods:
 * - NM_MIN_MAX: rescales values to the range [0, 1]
 * - NM_STANDARDIZE: applies z-score normalization (zero mean, unit variance)
 * - NM_ROBUST: applies robust scaling using median and interquartile range
 * (IQR)
 *
 * @param array  Input/output array to normalize in-place.
 * @param method Normalization strategy to apply.
 *
 * @warning For NM_ROBUST, division by a near-zero interquartile range may lead
 * to numerical instability if the input is nearly constant.
 *
 * See unit tests: @ref test_normalized.cpp
 */
void normalize(Array &array, NormalizationMethod method);

/**
 * @brief Returns a normalized copy of an array.
 *
 * This is a convenience wrapper around normalize().
 *
 * @param  array  Input array to normalize.
 * @param  method Normalization strategy to apply.
 *
 * @return        A new Array containing the normalized values.
 *
 * See unit tests: @ref test_normalized.cpp
 */
Array normalized(const Array &array, NormalizationMethod method);

/**
 * @brief Compute the population variance of a 2D array.
 *
 * @param  array  Input 2D array, accessed via array(i, j).
 * @param  p_mean Pointer to a precomputed mean (avoids recomputation); if null,
 *                the mean is computed internally.
 * @return        Population variance.
 *
 * See unit tests: @ref test_variance.cpp
 */
float variance(const Array &array, float *p_mean = nullptr);

} // namespace hmap