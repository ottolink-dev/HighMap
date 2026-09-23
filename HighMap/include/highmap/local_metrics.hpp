/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file local_metrics.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Computes the downstream path length to the domain outlet.
 *
 * For each cell, follows the steepest-descending flow path and returns the
 * number of downstream cells required to leave the domain. Cells whose flow
 * terminates in an internal sink have a path length of 0.
 *
 * **Example**
 * @include ex_path_length_to_sink.cpp
 *
 * **Result**
 * @image html ex_path_length_to_sink.png
 */
Array path_length_to_outlet(const Array &z);

/**
 * @brief Computes the downstream path length to the nearest sink.
 *
 * For each cell, follows the steepest-descending flow path and returns the
 * number of downstream cells required to reach a sink. Sink cells have a path
 * length of 0.
 *
 * **Example**
 * @include ex_path_length_to_sink.cpp
 *
 * **Result**
 * @image html ex_path_length_to_sink.png
 */
Array path_length_to_sink(const Array &z);

/**
 * @brief Return the local maxima based on a maximum filter with a square
 * kernel.
 *
 * This function identifies the local maxima in the input array using a maximum
 * filter with a square kernel. The local maxima are determined based on the
 * footprint radius specified by `ir`. The result is an array where each value
 * represents the local maximum within the defined kernel size.
 *
 * @param  array Input array from which local maxima are to be extracted.
 * @param  ir    Square kernel footprint radius. The size of the kernel used to
 *               determine the local maxima.
 * @return       Array Resulting array containing the local maxima.
 *
 * **Example**
 * @include ex_local_max.cpp
 *
 * **Result**
 * @image html ex_local_max.png
 *
 * @see          {@link local_max_disk}, {@link local_min}
 */
Array local_max(const Array &array, int ir);

/**
 * @brief Computes the local median deviation of a 2D array.
 *
 * This function calculates the absolute difference between the local mean and a
 * pseudo-local median of each element in the input array. The local
 * neighborhood is defined by a square window with radius `ir`.
 *
 * @param  array The input 2D array of values (e.g., a heightmap).
 * @param  ir    The radius of the square neighborhood window used for computing
 *               local statistics. The window size is (2*ir + 1) x (2*ir + 1).
 *
 * @return       A new array of the same size as `array`, where each element
 *               represents the absolute deviation between the local mean and
 *               pseudo-local median in its neighborhood.
 *
 * @note The function uses a pseudo-median approximation. For exact median
 * computation, replace the `median_pseudo()` call with a proper median filter
 * implementation.
 *
 * **Example**
 * @include ex_local_median_deviation.cpp
 *
 * **Result**
 * @image html ex_local_median_deviation.png
 *
 * @see          local_mean(), median_pseudo()
 */
Array local_median_deviation(const Array &array, int ir);

/**
 * @brief Return the local minima based on a maximum filter with a square
 * kernel.
 *
 * This function identifies the local minima in the input array using a maximum
 * filter with a square kernel. The local minima are determined based on the
 * footprint radius specified by `ir`. The result is an array where each value
 * represents the local minimum within the defined kernel size.
 *
 * @param  array Input array from which local minima are to be extracted.
 * @param  ir    Square kernel footprint radius. The size of the kernel used to
 *               determine the local minima.
 * @return       Array Resulting array containing the local minima.
 *
 * **Example**
 * @include ex_local_max.cpp
 *
 * **Result**
 * @image html ex_local_max.png
 *
 * @see          {@link local_min}
 */
Array local_min(const Array &array, int ir);

/**
 * @brief Return the local mean based on a mean filter with a square kernel.
 *
 * This function calculates the local mean of the input array using a mean
 * filter with a square kernel. The local mean is determined by averaging values
 * within a square neighborhood defined by the footprint radius `ir`. The result
 * is an array where each value represents the mean of the surrounding values
 * within the kernel size.
 *
 * @param  array Input array from which the local mean is to be calculated.
 * @param  ir    Square kernel footprint radius. The size of the kernel used to
 *               compute the local mean.
 * @return       Array Resulting array containing the local means.
 *
 * **Example**
 * @include ex_local_mean.cpp
 *
 * **Result**
 * @image html ex_local_mean0.png
 * @image html ex_local_mean1.png
 *
 * @see          {@link local_max}, {@link local_min}
 */
Array local_mean(const Array &array, int ir);

/**
 * @brief Calculates the relative elevation within a specified radius, helping
 * to identify local highs and lows.
 *
 * Relative elevation analysis determines how high or low a point is relative to
 * its surroundings, which can be critical in hydrological modeling and
 * geomorphology.
 *
 * **Usage**:
 * - Use this function to detect local depressions or peaks in the terrain,
 * which could indicate potential water collection points or hilltops.
 * - Useful in flood risk assessment and landscape classification.
 *
 * @param  array The input array representing the terrain elevation data.
 * @param  ir    The radius (in pixels) within which to calculate the relative
 *               elevation for each point.
 * @return       Array An output array containing the relative elevation values,
 *               normalized between 0 and 1.
 *
 * **Example**
 * @include ex_relative_elevation.cpp
 *
 * **Result**
 * @image html ex_relative_elevation.png
 */
Array relative_elevation(const Array &array, int ir);

/**
 * @brief Computes the ruggedness of each element in the input array.
 *
 * The ruggedness is calculated as the square root of the sum of squared
 * differences between each element and its neighbors within a specified radius.
 *
 * @param  array The input array for which ruggedness is to be computed.
 * @param  ir    The radius within which neighbors are considered for ruggedness
 *               calculation.
 * @return       An array containing the ruggedness values for each element in
 *               the input array.
 *
 * @note https://xdem.readthedocs.io/en/latest/terrain.html
 *
 * **Example**
 * @include ex_ruggedness.cpp
 *
 * **Result**
 * @image html ex_ruggedness0.png
 * @image html ex_ruggedness1.png
 */
Array ruggedness(const Array &array, int ir);

/**
 * @brief Estimates the rugosity of a surface by analyzing the skewness of the
 * elevation data, which reflects surface roughness.
 *
 * Rugosity is a measure of terrain roughness, often used in ecological studies
 * and habitat mapping. Higher rugosity values indicate more rugged terrain,
 * which can affect species distribution and water flow.
 *
 * @param  z      The input array representing the heightmap data (elevation
 *                values).
 * @param  ir     The radius of the square kernel used for calculations,
 *                determining the scale of the analysis.
 * @param  convex Return the convex rugosity if true, and the concave ones if
 *                not.
 * @return        Array An output array containing the rugosity estimates, where
 *                higher values indicate rougher terrain.
 *
 * **Example**
 * @include ex_rugosity.cpp
 *
 * **Result**
 * @image html ex_rugosity0.png
 * @image html ex_rugosity1.png
 */
Array rugosity(const Array &z, int ir, bool convex = true);

/**
 * @brief Measures the valley width by calculating the distance from each point
 * in a concave region to the frontier of that region.
 *
 * Valley width is a geomorphological metric that represents the distance across
 * valleys or concave regions in a terrain. This measurement is particularly
 * useful in hydrological modeling and landscape analysis, where valley
 * dimensions are important.
 *
 * **Usage**
 * - Use this function to assess the width of valleys within a terrain, which
 * can be important for understanding water flow, sediment transport, and
 * landform development.
 * - Applicable in studies focusing on river valleys, canyon analysis, and
 * erosion patterns.
 *
 * @param  z            The input array representing the heightmap data
 *                      (elevation values).
 * @param  ir           The radius used for pre-filtering, controlling the scale
 *                      of analysis (in pixels). The default value is 0, meaning
 *                      no pre-filtering is applied.
 * @param  ridge_select If enabled, selects ridges instead of valleys.
 * @return              Array An output array containing valley width values,
 *                      representing the distance to the edge of the concave
 *                      region for each point.
 *
 * **Example**
 * @include ex_valley_width.cpp
 *
 * **Result**
 * @image html ex_valley_width0.png
 * @image html ex_valley_width1.png
 * @image html ex_valley_width2.png
 */
Array valley_width(const Array &z, int ir = 0, bool ridge_select = false);

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Compute the average local aspect variance of an array.
 *
 * @param  array Input array.
 * @param  ir    Radius of the local window.
 * @return       Array containing the averaged local variance.
 *
 * **Example**
 * @include ex_local_metrics.cpp
 *
 * **Result**
 * @image html ex_local_metrics.png
 */
Array local_aspect_variance(const Array &array, int ir);

/**
 * @brief Compute the local maximum using a disk kernel.
 *
 * For each cell, returns the maximum value within a neighborhood of radius @p
 * ir defined by a disk-shaped kernel.
 *
 * @param  array Input array.
 * @param  ir    Radius of the disk kernel.
 * @return       Array of local maximum values.
 *
 *  **Example**
 * @include ex_local_metrics.cpp
 *
 * **Result**
 * @image html ex_local_metrics.png
 */
Array local_max(const Array &array, int ir);

/**
 * @brief Compute the local maximum using an octagonal kernel approximation
 * via 4 separable 1D passes (horizontal, vertical, and two diagonals).
 *
 * @param  array Input array.
 * @param  ir    Radius of the octagonal neighborhood footprint.
 * @return       Array of local maximum values.
 *
 * **Example**
 * @include ex_local_max_kernels.cpp
 *
 * **Result**
 * @image html ex_local_max_kernels.png
 */
Array local_max_octagon(const Array &array, int ir);

/**
 * @brief Compute the local maximum using a square kernel with a 2-pass
 * separable implementation.
 *
 * @param  array Input array.
 * @param  ir    Radius of the square kernel footprint.
 * @return       Array of local maximum values.
 *
 * **Example**
 * @include ex_local_max_kernels.cpp
 *
 * **Result**
 * @image html ex_local_max_kernels.png
 */
Array local_max_square(const Array &array, int ir);

/*! @brief See hmap::local_median_deviation */
Array local_median_deviation(const Array &array, int ir);

/**
 * @brief Compute the local minimum using a disk kernel.
 *
 * For each cell, returns the minimum value within a neighborhood of radius @p
 * ir defined by a disk-shaped kernel.
 *
 * @param  array Input array.
 * @param  ir    Radius of the disk kernel.
 * @return       Array of local minimum values.
 *
 *  **Example**
 * @include ex_local_metrics.cpp
 *
 * **Result**
 * @image html ex_local_metrics.png
 */
Array local_min(const Array &array, int ir);

/**
 * @brief Compute the local minimum using an octagonal kernel approximation
 * via 4 separable 1D passes (horizontal, vertical, and two diagonals).
 *
 * @param  array Input array.
 * @param  ir    Radius of the octagonal neighborhood footprint.
 * @return       Array of local minimum values.
 *
 * **Example**
 * @include ex_local_max_kernels.cpp
 *
 * **Result**
 * @image html ex_local_max_kernels.png
 */
Array local_min_octagon(const Array &array, int ir);

/**
 * @brief Compute the local minimum using a square kernel with a 2-pass
 * separable implementation.
 *
 * @param  array Input array.
 * @param  ir    Radius of the square kernel footprint.
 * @return       Array of local minimum values.
 *
 * **Example**
 * @include ex_local_max_kernels.cpp
 *
 * **Result**
 * @image html ex_local_max_kernels.png
 */
Array local_min_square(const Array &array, int ir);

/**
 * @brief Compute the local relief of an array.
 *
 * The local relief is defined as the difference between the maximum and minimum
 * values within a neighborhood of radius @p ir around each element. This
 * provides a measure of local variation in the array (e.g., terrain
 * ruggedness).
 *
 * @param  array Input array representing scalar values (e.g., elevation map).
 * @param  ir    Radius of the neighborhood (in pixels) used to compute local
 *               extrema.
 *
 * @return       Array An array of the same size as @p array, where each element
 *               contains the difference between the local maximum and minimum
 *               within the specified neighborhood.
 *
 * **Example**
 * @include ex_local_relief.cpp
 *
 * **Result**
 * @image html ex_local_relief.png
 */
Array local_relief(const Array &array, int ir);

/**
 * @brief Compute the average local variance of an array.
 *
 * Local variance is defined as the standard deviation within a window of radius
 * @p ir centered on each cell. The function computes this value for every cell,
 * then returns the average over the entire domain.
 *
 * @param  array Input array.
 * @param  ir    Radius of the local window.
 * @return       Array containing the averaged local variance.
 *
 * **Example**
 * @include ex_local_variance.cpp
 *
 * **Result**
 * @image html ex_local_variance.png
 */
Array local_variance(const Array &array, int ir);

/*! @brief See hmap::local_mean */
Array local_mean(const Array &array, int ir);

Array local_skewness(const Array &array, int ir);

/**
 * @brief Compute the local z-score of an array.
 *
 * For each cell, computes the normalized difference between its value and the
 * mean of a neighborhood of radius @p ir.
 *
 * @param  array Input array.
 * @param  ir    Radius of the local window.
 * @return       Array of local z-score values.
 *
 * **Example**
 * @include ex_local_z_score.cpp
 *
 * **Result**
 * @image html ex_local_z_score.png
 */
Array local_z_score(const Array &array, int ir);

/**
 * @brief Compute the Topographic Position Index (TPI).
 *
 * TPI measures the elevation difference between each cell and the mean
 * elevation within a neighborhood of radius @p ir.
 *
 * @param  array Input elevation array.
 * @param  ir    Radius of the neighborhood.
 * @return       Array of TPI values.
 *
 * **Example**
 * @include ex_topographic_position_index.cpp
 *
 * **Result**
 * @image html ex_topographic_position_index.png
 */
Array topographic_position_index(const Array &array, int ir);

/*! @brief See hmap::relative_elevation */
Array relative_elevation(const Array &array, int ir);

/*! @brief See hmap::relative_elevation */
Array relative_elevation_square_kernel(const Array &array, int ir);

/*! @brief See hmap::ruggedness */
Array ruggedness(const Array &array, int ir);

/*! @brief See hmap::rugosity */
Array rugosity(const Array &z, int ir, bool convex = true);

/*! @brief See hmap::valley_width */
Array valley_width(const Array &z, int ir = 0, bool ridge_select = false);

} // namespace hmap::gpu

// ==========================================================================
//  Wrapper
// ==========================================================================

namespace hmap::gpu
{

/**
 * \brief Types of local metrics.
 */
// clang-format off
enum LocalMetrics : int
{
	LM_LOCAL_ASPECT_VARIANCE,      ///< Variance of local slope directions.
	LM_LOCAL_MAX,                  ///< Local maximum value.
	LM_LOCAL_MEDIAN_DEVIATION,     ///< Deviation from local median.
	LM_LOCAL_MIN,                  ///< Local minimum value.
	LM_LOCAL_RELIEF,               ///< Diff. between local max & min.
	LM_LOCAL_VARIANCE,             ///< Variance of local values.
	LM_LOCAL_MEAN,                 ///< Mean of local values.
	LM_LOCAL_SKEWNESS,             ///< Skewness.
	LM_LOCAL_Z_SCORE,              ///< Standardized value.
	LM_TOPOGRAPHIC_POSITION_INDEX, ///< Topographic index.
	LM_RELATIVE_ELEVATION,         ///< Normalized elevation.
	LM_RUGGEDNESS,                 ///< Measure of terrain roughness.
	LM_RUGOSITY_CONCAVE,           ///< Roughness of concave features.
	LM_RUGOSITY_CONVEX,            ///< Roughness of convex features.
};
// clang-format on

/**
 * @brief Compute a local metric on the input array using a square neighborhood.
 *
 * Applies the selected local statistic or terrain descriptor (e.g. mean,
 * variance, ruggedness) over a neighborhood of radius \p ir.
 *
 * @param  array  Input array to analyze.
 * @param  ir     Radius of the square neighborhood.
 * @param  metric Local metric to compute.
 * @return        Array Output array containing the computed metric.
 *
 * **Example**
 * @include ex_local_metrics.cpp
 *
 * **Result**
 * @image html ex_local_metrics.png
 *
 * See unit tests: @ref test_local_metrics.cpp
 */
Array local_metrics(const Array &array, int ir, LocalMetrics metric);

} // namespace hmap::gpu