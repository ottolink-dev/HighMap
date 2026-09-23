/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file morphology.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for image morphology functions and operations.
 *
 * This header file provides declarations for functions and utilities related to
 * image morphology, including operations such as dilation, erosion, opening,
 * closing, etc. These functions are commonly used in image processing to
 * enhance or extract features based on the shape and structure of objects
 * within an image.
 *
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/kernels.hpp"
#include "highmap/local_metrics.hpp"

namespace hmap
{

/**
 * @brief Enumeration for different types of distance transforms.
 */
enum DistanceTransformType : int
{
  DT_EXACT,     ///< Exact distance transform.
  DT_APPROX,    ///< Approximate distance transform.
  DT_MANHATTAN, ///< Manhattan distance transform.
  DT_JFA,       ///< Approximate (JFA) distance transform.
};

/**
 * @brief Remove connected components smaller than a given size threshold.
 *
 * Detects connected regions in the input array and removes those whose area is
 * smaller than \p threshold_size.
 *
 * @param  array            Input array.
 * @param  threshold_size   Minimum component size to keep.
 * @param  background_value Value used to identify background.
 * @param  fill_value       Value used to replace removed components.
 * @return                  Array with small components removed.
 *
 * **Example**
 * @include ex_area_remove_fill.cpp
 *
 * **Result**
 * @image html ex_area_remove_fill.png
 *
 * See unit tests: @ref test_area_remove.cpp
 */
Array area_remove(const Array &array,
                  float        threshold_size,
                  float        background_value = 0.f,
                  float        fill_value = 0.f);

/**
 * @brief Apply a border algorithm to the input array using a square structure.
 *
 * @param  array Input array to which the border algorithm is applied.
 * @param  ir    Radius of the square kernel used for the border algorithm.
 * @return       Array Output array with the border applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array border(const Array &array, int ir);

/**
 * @brief Apply a closing algorithm to the input array using a square structure.
 *
 * @param  array Input array to which the closing algorithm is applied.
 * @param  ir    Radius of the square kernel used for the closing algorithm.
 * @return       Array Output array with the closing applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array closing(const Array &array, int ir);

/**
 * @brief Apply closing by reconstruction to the input array.
 *
 * Performs a dilation followed by a reconstruction by erosion, filling small
 * dark regions while preserving the shape of larger structures.
 *
 * @param  array        Input array.
 * @param  ir           Radius of the structuring element.
 * @param  k_smooth_max Smoothness factor (0 for hard max).
 * @return              Array Result after closing by reconstruction.
 */
Array closing_by_reconstruction(const Array &array,
                                int          ir,
                                float        k_smooth_max = 0.f);

/**
 * @brief Smooth contour boundaries of segmented regions.
 *
 * Reduces staircase artifacts on component borders by smoothing the transition
 * zone around edges.
 *
 * @param  array            Input labeled or binary array.
 * @param  ir               Radius of the smoothing kernel.
 * @param  transition_ratio Controls edge blending width.
 * @return                  Array with smoothed contours.
 *
 * **Example**
 * @include ex_contour_smoothing.cpp
 *
 * **Result**
 * @image html ex_contour_smoothing.png
 */
Array contour_smoothing(const Array &array,
                        int          ir,
                        float        transition_ratio = 0.1f);

/**
 * @brief Apply a dilation algorithm to the input array using a square
 * structure.
 *
 * @param  array Input array to which the dilation algorithm is applied.
 * @param  ir    Radius of the square kernel used for the dilation algorithm.
 * @return       Array Output array with the dilation applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array dilation(const Array &array, int ir);

/**
 * @brief Expand non-zero regions of an array by morphological dilation, while
 * preserving the original non-zero values.
 *
 * This function performs a local-maximum (dilation) operation on the input
 * array with a given radius @p ir, then copies back the original non-zero
 * values so they remain unchanged. Only the surrounding background is filled
 * with expanded values, effectively “growing” the borders of non-zero regions.
 *
 * @param  array Input 2D array where non-zero cells represent regions of
 *               interest and zero cells represent background.
 * @param  ir    Radius of the dilation (in grid cells). Defines how far the
 *               non-zero values can expand into the background.
 *
 * @return       A new 2D array where the borders of non-zero regions have been
 *               expanded outward by up to @p ir cells, while the original
 *               non-zero values are kept intact.
 */
Array dilation_expand_border_only(const Array &array, int ir);

Array dilation_expand_min_value_border_only(const Array &array);

/**
 * @brief Return the Euclidean distance transform.
 *
 * Exact transform based on Meijster et al. algorithm @cite Meijster2000.
 *
 * @param  array                   Input array to be transformed, will be
 *                                 converted into binary: 1 wherever input is
 *                                 greater than 0, 0 elsewhere.
 * @param  return_squared_distance Whether the distance returned is squared or
 *                                 not.
 * @return                         Array Reference to the output array.
 *
 * **Example**
 * @include ex_distance_transform.cpp
 *
 * **Result**
 * @image html ex_distance_transform0.png
 * @image html ex_distance_transform1.png
 * @image html ex_distance_transform2.png
 * @image html ex_distance_transform3.png
 * @image html ex_distance_transform4.png
 *
 * See unit tests: @ref test_distance_transform.cpp
 */
Array distance_transform(const Array &array,
                         bool         return_squared_distance = false);

/**
 * @brief Calculates an approximate distance transform of the input array.
 *
 * @param  array                   Input array to calculate the distance
 *                                 transform for.
 * @param  return_squared_distance Optional parameter to return squared
 *                                 distances. Defaults to false.
 * @return                         Array containing the approximate distance
 *                                 transform of the input array.
 *
 * **Example**
 * @include ex_distance_transform.cpp
 *
 * **Result**
 * @image html ex_distance_transform0.png
 * @image html ex_distance_transform1.png
 * @image html ex_distance_transform2.png
 * @image html ex_distance_transform3.png
 * @image html ex_distance_transform4.png
 *
 * See unit tests: @ref test_distance_transform.cpp
 */
Array distance_transform_approx(const Array &array,
                                bool         return_squared_distance = false);

/**
 * @brief Calculates the Manhattan distance transform of an array.
 *
 * @param  array                   Input array.
 * @param  return_squared_distance If true, returns the squared Manhattan
 *                                 distance instead of the actual distance.
 * @return                         Array containing the Manhattan distance
 *                                 transform of the input array.
 *
 * **Example**
 * @include ex_distance_transform.cpp
 *
 * **Result**
 * @image html ex_distance_transform0.png
 * @image html ex_distance_transform1.png
 * @image html ex_distance_transform2.png
 * @image html ex_distance_transform3.png
 * @image html ex_distance_transform4.png
 *
 * See unit tests: @ref test_distance_transform.cpp
 */
Array distance_transform_manhattan(const Array &array,
                                   bool return_squared_distance = false);

/**
 * @brief Return the Euclidean distance transform.
 *
 * Exact transform based on Meijster et al. algorithm @cite Meijster2000.
 *
 * @param  array                   Input array to be transformed, will be
 *                                 converted into binary: 1 wherever input is
 *                                 greater than 0, 0 elsewhere.
 * @param  return_squared_distance Whether the distance returned is squared or
 *                                 not.
 * @return                         Array Reference to the output array.
 *
 * See unit tests: @ref test_distance_transform.cpp
 */
Array distance_transform_with_closest(const Array     &array,
                                      Mat<glm::ivec2> &closest,
                                      bool return_squared_distance = false);

/**
 * @brief Apply an erosion algorithm to the input array using a square
 * structure.
 *
 * @param  array Input array to which the erosion algorithm is applied.
 * @param  ir    Radius of the square kernel used for the erosion algorithm.
 * @return       Array Output array with the erosion applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 *
 * See unit tests: @ref test_erosion.cpp
 */
Array erosion(const Array &array, int ir);

/**
 * @brief Apply a flood fill algorithm to the input array.
 *
 * @param array            Input array to be filled.
 * @param i                Seed point row index.
 * @param j                Seed point column index.
 * @param fill_value       Filling value.
 * @param background_value Background value.
 *
 * **Example**
 * @include ex_flood_fill.cpp
 *
 * **Result**
 * @image html ex_flood_fill.png
 */
void flood_fill(Array &array,
                int    i,
                int    j,
                float  fill_value = 1.f,
                float  background_value = 0.f);

/**
 * @brief Apply a morphological black hat algorithm to the input array using a
 * square structure.
 *
 * @param  array Input array to which the black hat algorithm is applied.
 * @param  ir    Radius of the square kernel used for the black hat algorithm.
 * @return       Array Output array with the black hat applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array morphological_black_hat(const Array &array, int ir);

/**
 * @brief Apply a morphological gradient algorithm to the input array using a
 * square structure.
 *
 * @param  array Input array to which the gradient algorithm is applied.
 * @param  ir    Radius of the square kernel used for the gradient algorithm.
 * @return       Array Output array with the gradient applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array morphological_gradient(const Array &array, int ir);

/**
 * @brief Apply a morphological Laplacian operator to the input array using a
 * square structure.
 *
 * Compared to the classical Laplacian, this operator is more robust to noise
 * and better preserves sharp features, as it relies on morphological extrema
 * rather than linear differences.
 *
 * @param  array Input array to which the morphological Laplacian is applied.
 * @param  ir    Radius of the square kernel used for the underlying dilation
 *               and erosion operations.
 * @return       Array Output array containing the morphological Laplacian.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array morphological_laplacian(const Array &array, int ir);

/**
 * @brief Apply a morphological top hat algorithm to the input array using a
 * square structure.
 *
 * @param  array Input array to which the top hat algorithm is applied.
 * @param  ir    Radius of the square kernel used for the top hat algorithm.
 * @return       Array Output array with the top hat applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array morphological_top_hat(const Array &array, int ir);

/**
 * @brief Apply an opening algorithm to the input array using a square
 * structure.
 *
 * @param  array Input array to which the opening algorithm is applied.
 * @param  ir    Radius of the square kernel used for the opening algorithm.
 * @return       Array Output array with the opening applied.
 *
 * **Example**
 * @include ex_morphology_base.cpp
 *
 * **Result**
 * @image html ex_morphology_base.png
 */
Array opening(const Array &array, int ir);

/**
 * @brief Apply opening by reconstruction to the input array.
 *
 * Performs an erosion followed by a reconstruction by dilation, removing small
 * bright features while preserving the shape of larger structures.
 *
 * @param  array        Input array.
 * @param  ir           Radius of the structuring element.
 * @param  k_smooth_min Smoothness factor (0 for hard min).
 * @return              Array Result after opening by reconstruction.
 */
Array opening_by_reconstruction(const Array &array,
                                int          ir,
                                float        k_smooth_min = 0.f);

/**
 * @brief Perform morphological reconstruction by dilation.
 *
 * Iteratively dilates the marker image while constraining values to be less
 * than or equal to the mask. The marker must be element-wise <= mask.
 * Optionally uses a smooth minimum controlled by \p k_smooth_min.
 *
 * @param  marker       Initial image (seed), must satisfy marker <= mask.
 * @param  mask         Constraint image (upper bound).
 * @param  ir           Radius of the structuring element.
 * @param  k_smooth_min Smoothness factor (0 for hard min).
 * @return              Array Reconstructed image.
 *
 * **Example**
 * @include ex_reconstruction_by_dilation.cpp
 *
 * **Result**
 * @image html ex_reconstruction_by_dilation.png
 */
Array reconstruction_by_dilation(const Array &marker,
                                 const Array &mask,
                                 int          ir,
                                 float        k_smooth_min = 0.f);

/**
 * @brief Perform morphological reconstruction by erosion.
 *
 * Iteratively erodes the marker image while constraining values to be greater
 * than or equal to the mask. The marker must be element-wise >= mask.
 * Optionally uses a smooth maximum controlled by \p k_smooth_max.
 *
 * @param  marker       Initial image (seed), must satisfy marker >= mask.
 * @param  mask         Constraint image (lower bound).
 * @param  ir           Radius of the structuring element.
 * @param  k_smooth_max Smoothness factor (0 for hard max).
 * @return              Array Reconstructed image.
 *
 * **Example**
 * @include ex_reconstruction_by_dilation.cpp
 *
 * **Result**
 * @image html ex_reconstruction_by_dilation.png
 */
Array reconstruction_by_erosion(const Array &marker,
                                const Array &mask,
                                int          ir,
                                float        k_smooth_max = 0.f);

/**
 * @brief Computes the relative distance of each non-zero cell in a binary array
 * from the skeleton and border.
 *
 * This function calculates a relative distance measure for each non-zero cell
 * in the input array. The measure is defined as the ratio of the cell's
 * distance to the nearest border and the combined distances to the nearest
 * skeleton and border cells. The skeleton is computed using the Zhang-Suen
 * skeletonization algorithm.
 *
 * @param  array           The input binary array for which the relative
 *                         distance map is to be calculated. Non-zero values are
 *                         considered for processing.
 * @param  ir_search       The search radius for finding the nearest skeleton
 *                         and border cells.
 * @param  zero_at_borders If true, the borders of the skeletonized image will
 *                         be set to zero.
 * @param  ir_erosion      The erosion radius applied to the skeleton.
 * @return                 An array representing the relative distance map,
 *                         where each cell has a value between 0 and 1. A value
 *                         closer to 1 indicates proximity to the skeleton,
 *                         while a value closer to 0 indicates proximity to the
 *                         border.
 *
 * @note The skeleton is computed using the Zhang-Suen skeletonization
 * algorithm.
 *
 * **Example**
 * @include ex_signed_curvature_from_distance.cpp
 *
 * **Result**
 * @image html ex_signed_curvature_from_distance0.png
 * @image html ex_signed_curvature_from_distance1.png
 * @image html ex_signed_curvature_from_distance2.png
 * @image html ex_signed_curvature_from_distance3.png
 */
Array relative_distance_from_skeleton(const Array &array,
                                      int          ir_search,
                                      bool         zero_at_borders = true,
                                      int          ir_erosion = 1);

Array relative_distance_from_skeleton(const Array &array,
                                      const Array &skeleton,
                                      int          ir_search,
                                      int ir_erosion = 1); ///< @overload

/**
 * @brief Computes a skeletonized version of an array.
 *
 * Extracts or refines the skeleton structure of the input array using a
 * reference skeleton mask. Optionally forces border values to zero to avoid
 * edge artifacts.
 *
 * @param  array           Input array to skeletonize.
 * @param  skeleton        Reference skeleton mask or initialization array.
 * @param  zero_at_borders If true, clears skeleton values at borders.
 * @return                 Skeletonized array.
 *
 * **Example**
 * @include ex_skeleton.cpp
 *
 * **Result**
 * @image html ex_skeleton.png
 */
Array skeleton(const Array &array,
               const Array &skeleton,
               bool         zero_at_borders = true);

/**
 * @brief Removes endpoint-like pixels from a binary/weighted array.
 *
 * Iteratively filters out isolated or weakly connected pixels based on their
 * neighborhood connectivity. Pixels matching the background value or having
 * fewer than 2 non-background neighbors are removed.
 *
 * @param  array            Input array to process.
 * @param  iterations       Number of erosion-like iterations to apply.
 * @param  background_value Value considered as background (removed).
 * @return                  Filtered array after endpoint removal.
 *
 * **Example**
 * @include ex_remove_endpoints.cpp
 *
 * **Result**
 * @image html ex_remove_endpoints.png
 */
Array remove_endpoints(const Array &array,
                       int          iterations = 1,
                       float        background_value = 0.f);

/**
 * @brief Computes the signed curvature of the distance transform.
 *
 * Computes the Euclidean distance transform of a binary array, optionally
 * smooths it, then returns the divergence of the normalized gradient: κ = div(
 * ∇d / ||∇d|| ). Positive values indicate convex regions, negative values
 * concave regions.
 *
 * @param  array        Binary input (non-zero = foreground).
 * @param  prefilter_ir Optional smoothing radius (0 = no smoothing).
 * @return              Curvature field of the distance transform.
 *
 * **Example**
 * @include ex_signed_curvature_from_distance.cpp
 *
 * **Result**
 * @image html ex_signed_curvature_from_distance0.png
 * @image html ex_signed_curvature_from_distance1.png
 * @image html ex_signed_curvature_from_distance2.png
 * @image html ex_signed_curvature_from_distance3.png
 */
Array signed_curvature_from_distance(const Array &array, int prefilter_ir = 0);

/**
 * @brief Computes a signed distance transform using curvature sign.
 *
 * First computes the unsigned Euclidean distance transform, then assigns its
 * sign from the curvature of the smoothed distance field: sign = sign( div( ∇d
 * / ||∇d|| ) ).
 *
 * @param  array        Binary input (non-zero = foreground).
 * @param  prefilter_ir Optional smoothing radius for curvature evaluation.
 * @return              Signed distance transform.
 */
Array signed_distance_transform(const Array &array, int prefilter_ir = 0);

/**
 * @brief Computes the skeleton of a binary image using the Zhang-Suen
 * skeletonization algorithm.
 *
 * This function processes a binary input array to extract its skeleton by
 * iteratively thinning the image until no further changes occur. It optionally
 * sets the borders of the resulting skeletonized image to zero.
 *
 * @param  array           The input binary array to be skeletonized. Values
 *                         should typically be 0 or 1.
 * @param  zero_at_borders If true, the borders of the resulting array will be
 *                         set to zero.
 * @return                 The skeletonized version of the input array.
 *
 * @note This implementation is based on the algorithm described at
 * https://github.com/krishraghuram/Zhang-Suen-Skeletonization.
 *
 * **Example**
 * @include ex_skeleton.cpp
 *
 * **Result**
 * @image html ex_skeleton.png
 *
 * See unit tests: @ref test_skeleton.cpp
 */
Array skeleton(const Array &array, bool zero_at_borders = true);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::border */
Array border(const Array &array,
             int          ir,
             MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::closing */
Array closing(const Array &array,
              int          ir,
              MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::closing_by_reconstruction */
Array closing_by_reconstruction(const Array &array,
                                int          ir,
                                float        k_smooth_max = 0.f,
                                MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::contour_smoothing */
Array contour_smoothing(const Array &array,
                        int          ir,
                        float        transition_ratio = 0.1f);

/*! @brief See hmap::dilation */
Array dilation(const Array &array,
               int          ir,
               MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::dilation_expand_border_only */
Array dilation_expand_border_only(
    const Array &array,
    int          ir,
    MinMaxKernel kernel_type = MinMaxKernel::DISK);

/**
 * @brief Return the Euclidean distance transform.
 *
 * (Almost) exact transform based on the jump flooding algorithm.
 *
 * @param  array                   Input array to be transformed, will be
 *                                 converted into binary: 1 wherever input is
 *                                 greater than 0, 0 elsewhere.
 * @param  return_squared_distance Whether the distance returned is squared or
 *                                 not.
 * @return                         Array Reference to the output array.
 *
 * **Example**
 * @include ex_distance_transform.cpp
 *
 * **Result**
 * @image html ex_distance_transform0.png
 * @image html ex_distance_transform1.png
 * @image html ex_distance_transform2.png
 * @image html ex_distance_transform3.png
 * @image html ex_distance_transform4.png
 */
Array distance_transform_jfa(const Array &array,
                             bool         return_squared_distance = false);

/*! @brief See hmap::erosion */
Array erosion(const Array &array,
              int          ir,
              MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::morphological_black_hat */
Array morphological_black_hat(const Array &array,
                              int          ir,
                              MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::morphological_gradient */
Array morphological_gradient(const Array &array,
                             int          ir,
                             MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::morphological_laplacian */
Array morphological_laplacian(const Array &array,
                              int          ir,
                              MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::morphological_top_hat */
Array morphological_top_hat(const Array &array,
                            int          ir,
                            MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::opening */
Array opening(const Array &array,
              int          ir,
              MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::opening_by_reconstruction */
Array opening_by_reconstruction(const Array &array,
                                int          ir,
                                float        k_smooth_min = 0.f,
                                MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::reconstruction_by_dilation */
Array reconstruction_by_dilation(const Array &marker,
                                 const Array &mask,
                                 int          ir,
                                 float        k_smooth_min = 0.f,
                                 MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::reconstruction_by_erosion */
Array reconstruction_by_erosion(const Array &marker,
                                const Array &mask,
                                int          ir,
                                float        k_smooth_max = 0.f,
                                MinMaxKernel kernel_type = MinMaxKernel::DISK);

/*! @brief See hmap::relative_distance_from_skeleton */
Array relative_distance_from_skeleton(const Array &array,
                                      int          ir_search,
                                      bool         zero_at_borders = true,
                                      int          ir_erosion = 1);

/*! @brief See hmap::relative_distance_from_skeleton */
Array relative_distance_from_skeleton(const Array &array,
                                      const Array &skeleton,
                                      int          ir_search,
                                      int          ir_erosion = 1);

/*! @brief See hmap::signed_curvature_from_distance */
Array signed_curvature_from_distance(const Array &array, int prefilter_ir = 0);

/*! @brief See hmap::signed_distance_transform */
Array signed_distance_transform(const Array &array, int prefilter_ir = 0);

/*! @brief See hmap::skeleton */
Array skeleton(const Array &array, bool zero_at_borders = true);

} // namespace hmap::gpu

// ==========================================================================
//  Wrapper
// ==========================================================================

namespace hmap
{

/**
 * \brief Types of base morphology operators.
 */
// clang-format off
enum MorphologyOperation : int
{
	MO_BORDER,    //!< Extract border of structures.
	MO_CLOSING,   //!< Dilation followed by erosion (fills small holes).
	MO_DILATION,  //!< Expand structures (local maximum).
	MO_EROSION,   //!< Shrink structures (local minimum).
	MO_OPENING,   //!< Erosion followed by dilation (removes small objects).
	MO_BLACK_HAT, //!< Difference between closing and input (dark features).
	MO_TOP_HAT,   //!< Difference between input and opening (bright
	// features).
	MO_GRADIENT,  //!< Difference between dilation and erosion (edge
	// magnitude).
	MO_LAPLACIAN, //!< Second-order operator highlighting ridges and
	// valleys.
	MO_CLOSING_BY_RECONSTRUCTION,
	MO_OPENING_BY_RECONSTRUCTION,
};
// clang-format on

/**
 * @brief Apply a morphological operation to the input array using a square
 * kernel.
 *
 * Dispatches to the selected morphological operator (e.g. erosion, dilation,
 * opening, closing, gradient, etc.).
 *
 * @param  array     Input array to process.
 * @param  ir        Radius of the square structuring element.
 * @param  operation Morphological operation to apply.
 * @return           Array Resulting array after applying the operation.
 */
Array morphological_operators(const Array        &array,
                              int                 ir,
                              MorphologyOperation operation);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::morphological_operators */
Array morphological_operators(const Array        &array,
                              int                 ir,
                              MorphologyOperation operation);

} // namespace hmap::gpu