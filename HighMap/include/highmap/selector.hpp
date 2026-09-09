/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file op.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for region selection functions in 2D array data
 * representation.
 *
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/range.hpp"

namespace hmap
{

/**
 * @brief Perturb the contour of a binary mask using a displacement field, while
 * trying to preserve a single connected component without holes.
 *
 * This function takes a binary mask and a noise map of the same size, computes
 * the boundary of the mask, estimates a local outward normal for each boundary
 * pixel, and displaces that pixel by an amount proportional to the noise value.
 * The displacement is clamped to remain inside the image domain.
 *
 * @param  mask             Input binary mask as an Array (0 = background, >0 =
 *                          foreground). Must be at least 3×3 in size.
 * @param  noise            Noise field providing signed displacement factors,
 *                          typically in [-1, 1]. Must have the same dimensions
 *                          as @p mask.
 * @param  max_displacement Maximum displacement magnitude in pixels applied
 *                          along the estimated boundary normal for each pixel.
 * @param  ir               Neighborhood radius used for normal estimation and
 *                          contour reconnection.
 *
 * @return                  A new Array representing the perturbed binary mask,
 *                          with values 0 (background) and 1 (foreground),
 *                          guaranteed to be a single connectedcomponent without
 *                          holes.
 *
 * **Example**
 * @include ex_perturb_mask_contour.cpp
 *
 * **Result**
 * @image html ex_perturb_mask_contour.png
 */
Array perturb_mask_contour(const Array &mask,
                           const Array &noise,
                           float        max_displacement,
                           int          ir = 1);

/**
 * @brief Mask adjustement using a 'scanning' method.
 *
 * See https://www.shadertoy.com/view/stjSRR
 *
 * @param  array      Input array.
 * @param  contrast   Contrast.
 * @param  brightness Brightness.
 * @return            Ouput array.
 *
 * **Example**
 * @include ex_scan_mask.cpp
 *
 * **Result**
 * @image html ex_scan_mask.png
 */
Array scan_mask(const Array &array,
                float        contrast = 0.5f,
                float        brightness = 0.5f);

/**
 * @brief Return angle selection for a given angle and a tolerance half-width on
 * this value.
 * @param  array Input array.
 * @param  angle Selected angle (degree).
 * @param  sigma Selected angle tolerance (degree).
 * @param  ir    Prefilter radius.
 * @return       Output array.
 *
 * **Example**
 * @include ex_select_angle.cpp
 *
 * **Result**
 * @image html ex_select_angle.png
 */
Array select_angle(const Array &array, float angle, float sigma, int ir = 0);

/**
 * @brief Return blob detection using the Laplacian of Gaussian (LoG) approach.
 *
 * @param  array Input array.
 * @param  ir    Kernel radius.
 * @return       Array Output array.
 *
 * **Example**
 * @include ex_select_blob.cpp
 *
 * **Result**
 * @image html ex_select_blob0.png
 * @image html ex_select_blob1.png
 */
Array select_blob_log(const Array &array, int ir);

/**
 * @brief Return holes or bumps detection based on the mean curvature of the
 * heightmap.
 *
 * @param  array   Input array.
 * @param  ir      Kernel radius.
 * @param  concave Select 'holes' if set to true, and select 'bumps' if set to
 *                 false.
 * @return         Array Output array.
 *
 * **Example**
 * @include ex_select_cavities.cpp
 *
 * **Result**
 * @image html ex_select_cavities.png
 */
Array select_cavities(const Array &array, int ir, bool concave = true);

/**
 * @brief
 *
 * @param  array
 * @return       Array
 *
 * **Example**
 * @include ex_select_elevation_slope.cpp
 *
 * **Result**
 * @image html ex_select_elevation_slope0.png
 * @image html ex_select_elevation_slope1.png
 */
Array select_elevation_slope(const Array &array, float gradient_scale);
Array select_elevation_slope(const Array &array,
                             float        gradient_scale,
                             float        vmax); ///< @overload

/**
 * @brief Return an array with elements equal to 1 where input elements are
 * equal to `value`.
 *
 * @param  array Input array.
 * @param  value Criteria value.
 * @return       Array Output array.
 *
 * **Example**
 * @include ex_select.cpp
 *
 * **Result**
 * @image html ex_select0.png
 * @image html ex_select1.png
 * @image html ex_select2.png
 * @image html ex_select3.png
 * @image html ex_select4.png
 */
Array select_eq(const Array &array, float value);

/**
 * @brief Return an array weighted by the gap between the gradient angle and a
 * given angle.
 *
 * @param  array Input array.
 * @param  angle Reference angle (degree).
 * @return       Array Output array.
 */
Array select_gradient_angle(const Array &array, float angle);

/**
 * @brief Return an array filled with 1 where the gradient is larger than a
 * given value and 0 elsewhere.
 *
 * @param  array        Input array.
 * @param  talus_center Reference talus.
 * @return              Array Output array.
 */
Array select_gradient_binary(const Array &array, float talus_center);

/**
 * @brief Return an array weighted (exponantial decay) by the gradient norm of
 * the input array.
 *
 * @param  array        Input array.
 * @param  talus_center Reference talus.
 * @param  talus_sigma  Talus half-width.
 * @return              Array Output array.
 */
Array select_gradient_exp(const Array &array,
                          float        talus_center,
                          float        talus_sigma);

/**
 * @brief Return an array weighted (square inverse) by the gradient norm of the
 * input array.
 *
 * @param  array        Input array.
 * @param  talus_center Reference talus.
 * @param  talus_sigma  Talus half-width.
 * @return              Array Output array.
 */
Array select_gradient_inv(const Array &array,
                          float        talus_center,
                          float        talus_sigma);

/**
 * @brief Return an array with elements equal to 1 where input elements are
 * larger than `value`.
 *
 * @param  array Input array.
 * @param  value Criteria value.
 * @return       Array Output array.
 *
 * **Example**
 * @include ex_select.cpp
 *
 * **Result**
 * @image html ex_select0.png
 * @image html ex_select1.png
 * @image html ex_select2.png
 * @image html ex_select3.png
 * @image html ex_select4.png
 */
Array select_gt(const Array &array, float value);

/**
 * @brief Return an array with elements equal to 1 where input elements are
 * within the bounds provided.
 *
 * @param  array  Input array.
 * @param  value1 Lower bound.
 * @param  value2 Upper bound.
 * @return        Array Output array.
 *
 * **Example**
 * @include ex_select.cpp
 *
 * **Result**
 * @image html ex_select0.png
 * @image html ex_select1.png
 * @image html ex_select2.png
 * @image html ex_select3.png
 * @image html ex_select4.png
 */
Array select_interval(const Array &array, float value1, float value2);

/**
 * @brief Return an array with positive values if the slope is pointing to the
 * center (slope is inward), and negative values otherwise (slope is outward).
 * @param  array  Input array.
 * @param  center Reference center.
 * @param  bbox   Domain bounding box.
 * @return        Output array.
 *
 * **Example**
 * @include ex_select_inward_outward_slope.cpp
 *
 * **Result**
 * @image html ex_select_inward_outward_slope.png
 */
Array select_inward_outward_slope(const Array &array,
                                  glm::vec2    center = {0.5f, 0.5f},
                                  glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array with elements equal to 1 where input elements are
 * smaller than `value`.
 *
 * @param  array Input array.
 * @param  value Criteria value.
 * @return       Array Output array.
 *
 * **Example**
 * @include ex_select.cpp
 *
 * **Result**
 * @image html ex_select0.png
 * @image html ex_select1.png
 * @image html ex_select2.png
 * @image html ex_select3.png
 * @image html ex_select4.png
 */
Array select_lt(const Array &array, float value);

/**
 * @brief Selects the midrange values of the input array within a specified
 * range.
 *
 * @param  array The input array from which the midrange values are selected.
 * @param  gain  The gain factor to scale the selected midrange values.
 * @param  vmin  The minimum value for range.
 * @param  vmax  The maximum value for range.
 * @return       An array containing the midrange values.
 *
 * **Example**
 * @include ex_select_midrange.cpp
 *
 * **Result**
 * @image html ex_select_midrange.png
 */
Array select_midrange(const Array &array, float gain, float vmin, float vmax);

/**
 * @brief Selects and scales the midrange values of the input array.
 *
 * This function extracts the midrange portion of the input array and scales it
 * by the specified gain.
 *
 * @param  array The input array from which the midrange values are selected.
 * @param  gain  The gain factor to scale the selected midrange values.
 * @return       An array containing the midrange values scaled by the gain.
 *
 * **Example**
 * @include ex_select_midrange.cpp
 *
 * **Result**
 * @image html ex_select_midrange.png
 */
Array select_midrange(const Array &array, float gain);

/**
 * @brief Splits the input array into three bands (low, mid, and high) based on
 * given ratios and overlap.
 *
 * This function processes the input array and divides it into three value
 * bands: low, mid, and high. The bands are determined by the provided ratios
 * and overlap. Optionally, a range (vmin, vmax) can be applied to the values in
 * the array.
 *
 * @param array     The input array that contains the data to be split into
 *                  bands.
 * @param band_low  The output array for the low-value band.
 * @param band_mid  The output array for the mid-value band.
 * @param band_high The output array for the high-value band.
 * @param ratio1    The ratio that defines the split between the low and mid
 *                  bands.
 * @param ratio2    The ratio that defines the split between the mid and high
 *                  bands.
 * @param overlap   The amount of overlap between adjacent bands.
 * @param vmin      The minimum value for range filtering. Values below vmin
 *                  will be clamped. (Optional, only in the first function)
 * @param vmax      The maximum value for range filtering. Values above vmax
 *                  will be clamped. (Optional, only in the first function)
 *
 * **Example**
 * @include ex_select_multiband3.cpp
 *
 * **Result**
 * @image html ex_select_multiband3.png
 */
void select_multiband3(const Array &array,
                       Array       &band_low,
                       Array       &band_mid,
                       Array       &band_high,
                       float        ratio1,
                       float        ratio2,
                       float        overlap,
                       float        vmin,
                       float        vmax);

/**
 * @brief Splits the input array into three bands (low, mid, and high) based on
 * given ratios and overlap.
 *
 * This function processes the input array and divides it into three value
 * bands: low, mid, and high. The bands are determined by the provided ratios
 * and overlap.
 *
 * @param array     The input array that contains the data to be split into
 *                  bands.
 * @param band_low  The output array for the low-value band.
 * @param band_mid  The output array for the mid-value band.
 * @param band_high The output array for the high-value band.
 * @param ratio1    The ratio that defines the split between the low and mid
 *                  bands.
 * @param ratio2    The ratio that defines the split between the mid and high
 *                  bands.
 * @param overlap   The amount of overlap between adjacent bands.
 *
 * **Example**
 * @include ex_select_multiband3.cpp
 *
 * **Result**
 * @image html ex_select_multiband3.png
 */
void select_multiband3(const Array &array,
                       Array       &band_low,
                       Array       &band_mid,
                       Array       &band_high,
                       float        ratio1,
                       float        ratio2,
                       float        overlap);

/**
 * @brief Return an array filled with non-zero values where the input is in the
 * interval [value - sigma, value + sigma]. Output array values have a cubic
 * pulse distribution within this interval.
 *
 * @param  array Input array.
 * @param  value Central value.
 * @param  sigma Pulse width.
 * @return       Array Output array.
 *
 * **Example**
 * @include ex_select_pulse.cpp
 *
 * **Result**
 * @image html ex_select_pulse.png
 */
Array select_pulse(const Array &array, float value, float sigma);

/**
 * @brief Return an array filled with a criterion based on the occurence of a
 * river bed.
 *
 * @param  array          Input array.
 * @param  talus_ref      Reference talus used to localy define the
 *                        flow-partition.
 * @param  clipping_ratio Flow accumulation clipping ratio.
 * @return                Array Output array.
 *
 * **Example**
 * @include ex_select_rivers.cpp
 *
 * **Result**
 * @image html ex_select_rivers.png
 */
Array select_rivers(const Array &array, float talus_ref, float clipping_ratio);

/**
 * @brief Return an array filled with 1 at the blending transition between two
 * arrays, and 0 elsewhere.
 *
 * @param  array1      Input array 1.
 * @param  array2      Input array 2.
 * @param  array_blend Blending of array 1 and 2 to analyze.
 * @return             Array Resulting array
 *
 * **Example**
 * @include ex_select_transitions.cpp
 *
 * **Result**
 * @image html ex_select_transitions0.png
 * @image html ex_select_transitions1.png
 * @image html ex_select_transitions2.png
 */
Array select_transitions(const Array &array1,
                         const Array &array2,
                         const Array &array_blend);

Array select_valley(const Array &z, int ir, bool ridge_select = false);

/**
 * @brief Smooths a binary mask with an inward falloff while preserving its
 * 0-boundary.
 *
 * Computes a smooth falloff inside the non-zero region of the mask towards the
 * 0/1 frontier using distance transform and smoothstep. The 0-set of the output
 * remains identical to the input mask (no dilation beyond the boundary).
 *
 * @param  mask          Input binary mask ({0, 1}).
 * @param  feather_width Width of the feathering falloff in grid units.
 * @return               Array         Smoothed mask with smooth falloff of the
 *                       1s.
 */
Array smooth_mask_preserve_frontier(const Array &mask, float feather_width);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::select_cavities */
Array select_cavities(const Array &array, int ir, bool concave = true);

/**
 * @brief Computes a soil–flow selection map based on terrain gradient, river
 * mask, and smoothing parameters.
 *
 * @param  z                       Input heightmap.
 * @param  ir_gradient             Radius used for the morphological gradient.
 * @param  gradient_weight         Weight applied to the gradient magnitude.
 * @param  gradient_scaling_factor Scaling factor applied to the gradient;
 *                                 defaults to z.shape.x if <= 0.
 * @param  talus_ref               Reference talus value for river extraction;
 *                                 defaults to 1.f / z.shape.x if <= 0.
 * @param  clipping_ratio          Clipping ratio used in the river selection
 *                                 step.
 * @param  flow_gamma              Gamma correction exponent applied to the
 *                                 river mask; 1.0 disables correction.
 * @param  k_smooth                Smoothing parameter for the maximum blending.
 *
 * @return                         A 2D Array representing the final soil-flow
 *                                 selection mask.
 *
 * **Example**
 * @include ex_select_soil_weathered.cpp
 *
 * **Result**
 * @image html ex_select_soil_weathered.png
 */
Array select_soil_flow(const Array &z,
                       int          ir_gradient = 1,
                       float        gradient_weight = 1.f,
                       float        gradient_scaling_factor = 0.f,
                       float        flow_weight = 0.05f,
                       float        talus_ref = 0.f,
                       float        clipping_ratio = 50.f,
                       float        flow_gamma = 1.f,
                       float        k_smooth = 0.01f);

/**
 * @brief Computes a multi-scale soil/rock selector using curvature analysis.
 *
 * Applies progressively increasing smoothing radii in logarithmic steps
 * (handling ir_min = 0), computes the mean curvature at each scale, clamps it,
 * and accumulates the maximum response across all scales.
 *
 * @param  z                     Input heightmap or scalar field.
 * @param  ir_max                Maximum smoothing radius (log-scale upper
 *                               bound).
 * @param  ir_min                Minimum smoothing radius (0 allowed; replaced
 *                               by 1 for log scale).
 * @param  steps                 Number of logarithmic scales to evaluate.
 * @param  smaller_scales_weight Weight factor applied when descending scales.
 * @param  curvature_clamp_mode  Mode used to clamp curvature values.
 * @param  curvature_clamping    Curvature clamp threshold.
 *
 * @return                       An Array representing soil/rock selection
 *                               values across scales.
 *
 * **Example**
 * @include ex_select_soil_rocks.cpp
 *
 * **Result**
 * @image html ex_select_soil_rocks.png
 */

Array select_soil_rocks(
    const Array &z,
    int          ir_max = 64,
    int          ir_min = 0,
    int          steps = 4,
    float        smaller_scales_weight = 1.f,
    ClampMode    curvature_clamp_mode = ClampMode::POSITIVE_ONLY,
    float        curvature_clamping = 1.f);

/**
 * @brief Computes a soil weathering selection map based on curvature and
 * gradient analysis.
 *
 * This function generates a combined array representing the degree of soil
 * weathering using two main components:
 * - **Curvature contribution**: derived from smoothed height variations.
 * - **Gradient contribution**: derived from local slope magnitudes.
 *
 * @param z                       Input heightmap or scalar field.
 * @param ir_curvature            Radius (in pixels or samples) for curvature
 *                                smoothing. If zero, curvature smoothing is
 *                                skipped.
 * @param ir_gradient             Radius (in pixels or samples) for gradient
 *                                computation.
 * @param curvature_clamp_mode    Clamping mode applied to curvature values (see
 *                                ::ClampMode).
 * @param curvature_clamping      Maximum absolute value used for curvature
 *                                clamping.
 * @param curvature_weight        Weight applied to the curvature component in
 *                                the final mix.
 * @param gradient_weight         Weight applied to the gradient component in
 *                                the final mix.
 * @param gradient_scaling_factor Empirical scaling factor for both curvature
 *                                and gradient terms. If non-positive, defaults
 *                                to `z.shape.x`.
 *
 * **Example**
 * @include ex_select_soil_weathered.cpp
 *
 * **Result**
 * @image html ex_select_soil_weathered.png
 */
Array select_soil_weathered(
    const Array &z,
    int          ir_curvature = 0,
    int          ir_gradient = 4,
    ClampMode    curvature_clamp_mode = ClampMode::POSITIVE_ONLY,
    float        curvature_clamping = 10.f,
    float        curvature_weight = 1.f,
    float        gradient_weight = 1.f,
    float        gradient_scaling_factor = 0.f);

// if the gradient norm is known before-hand
Array select_soil_weathered(const Array &z,
                            const Array &gradient_norm, // in [0, 1]
                            int          ir_curvature,
                            ClampMode    curvature_clamp_mode,
                            float        curvature_clamping,
                            float        curvature_weight,
                            float        gradient_weight,
                            float gradient_scaling_factor); // for curvature

/*! @brief See hmap::select_valley */
Array select_valley(const Array &z, int ir, bool ridge_select = false);

} // namespace hmap::gpu
