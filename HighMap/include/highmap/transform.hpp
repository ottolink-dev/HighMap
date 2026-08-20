/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file transform.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for array transformation functions including rotation,
 * flipping, warping, and transposition.
 *
 * This header file provides declarations for various functions used to
 * transform arrays in different ways. It includes functionalities for rotating
 * arrays, flipping them horizontally or vertically, applying warping effects,
 * and transposing arrays. These transformations are useful in image processing,
 * data manipulation, and other applications requiring geometric modifications
 * of arrays.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include "macrologger.h"

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Flip the array horizontally (left/right).
 *
 * This function flips the input array along the vertical axis, resulting in a
 * left-to-right mirror image of the original array.
 *
 * @param array Input array to be flipped horizontally.
 *
 * **Example**
 * @include ex_flip_ud.cpp
 *
 * **Result**
 * @image html ex_flip_ud.png
 */
void flip_lr(Array &array);

/**
 * @brief Flip the array vertically (up/down).
 *
 * This function flips the input array along the horizontal axis, resulting in
 * an up-to-down mirror image of the original array.
 *
 * @param array Input array to be flipped vertically.
 *
 * **Example**
 * @include ex_flip_ud.cpp
 *
 * **Result**
 * @image html ex_flip_ud.png
 */
void flip_ud(Array &array);

/**
 * @brief Interpret the input array `dr` as a radial displacement and convert it
 * to a pair of displacements `dx` and `dy` in cartesian coordinates.
 * @param dr        Radial displacement.
 * @param dx        Displacent in x direction (output).
 * @param dy        Displacent in y direction (output).
 * @param smoothing Smoothing parameter to avoid discontinuity at the origin.
 * @param center    Origin center.
 * @param bbox      Domain bounding box.
 *
 * **Example**
 * @include ex_radial_displacement_to_xy.cpp
 *
 * **Result**
 * @image html ex_radial_displacement_to_xy.png
 * */
void radial_displacement_to_xy(const Array &dr,
                               Array       &dx,
                               Array       &dy,
                               float        smoothing = 1.f,
                               glm::vec2    center = {0.5f, 0.5f},
                               glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Rotate the array by 180 degrees.
 *
 * This function rotates the input array by 180 degrees in the counterclockwise
 * direction. The dimensions of the array will be adjusted accordingly.
 *
 * @param array Input array to be rotated by 180 degrees.
 */
void rot180(Array &array);

/**
 * @brief Rotate the array by 270 degrees.
 *
 * This function rotates the input array by 270 degrees in the counterclockwise
 * direction. The dimensions of the array will be adjusted accordingly.
 *
 * @param array Input array to be rotated by 270 degrees.
 */
void rot270(Array &array);

/**
 * @brief Rotate the array by 90 degrees.
 *
 * This function rotates the input array by 90 degrees in the counterclockwise
 * direction. The dimensions of the array will be adjusted accordingly.
 *
 * @param array Input array to be rotated by 90 degrees.
 *
 * **Example**
 * @include ex_rot90.cpp
 *
 * **Result**
 * @image html ex_rot90.png
 */
void rot90(Array &array);

/**
 * @brief Rotate the array by a specified angle.
 *
 * This function rotates the input array by a given angle in degrees. The
 * rotation can be performed with optional zero-padding, which fills in the
 * borders with zeros instead of using symmetric padding. This is particularly
 * useful when the array contains zero values at its borders.
 *
 * @param array        Input array to be rotated.
 * @param angle        Rotation angle in degrees.
 * @param zoom_in      If true, zoom in after rotation (default is true).
 * @param zero_padding If true, use zero-padding to fill the borders; otherwise,
 *                     use symmetry (default is false).
 *
 * **Example**
 * @include ex_rotate.cpp
 *
 * **Result**
 * @image html ex_rotate.png
 */
void rotate(Array &array,
            float  angle,
            bool   zoom_in = true,
            bool   zero_padding = false);

/**
 * @brief Rotates a scalar displacement field into directional X and Y
 * components.
 *
 * This function converts a scalar displacement map into two directional
 * components (`dx` and `dy`) according to a given rotation angle. It is
 * typically used to orient procedural displacements, such as those applied to
 * terrain heightmaps or noise-based distortion fields.
 *
 * @param delta Input scalar displacement field.
 * @param angle Rotation angle in degrees (counterclockwise).
 * @param dx    Output array receiving the X-axis displacement component.
 * @param dy    Output array receiving the Y-axis displacement component.
 * */
void rotate_displacement(const Array &delta, float angle, Array &dx, Array &dy);

/**
 * @brief Return the transposed array.
 *
 * This function returns a new array that is the transpose of the input array.
 * The transpose operation swaps the rows and columns of the array, effectively
 * flipping the array over its diagonal.
 *
 * @param  array Input array to be transposed.
 * @return       Array The transposed array.
 */
Array transpose(const Array &array);

/**
 * @brief Translates a 2D array by a specified amount along the x and y axes.
 *
 * This function shifts the contents of the input array by `dx` and `dy` units
 * along the x and y axes, respectively. It supports both periodic boundary
 * conditions, where the array wraps around, and non-periodic conditions, where
 * the shifted areas are filled with zeros.
 *
 * @param  array     The input 2D array to be translated. This array remains
 *                   unmodified.
 * @param  dx        The translation distance along the x-axis. Positive values
 *                   shift the array to the right.
 * @param  dy        The translation distance along the y-axis. Positive values
 *                   shift the array downward.
 * @param  periodic  If set to `true`, the translation is periodic, meaning that
 *                   elements that move out of one side of the array reappear on
 *                   the opposite side. If `false`, the areas exposed by the
 *                   translation are filled with zeros. The default is `false`.
 * @param  p_noise_x Optional pointer to a 2D array that contains x-direction
 *                   noise to be added to the translation. If provided, the
 *                   noise values are added to `dx` on a per-element basis.
 * @param  p_noise_y Optional pointer to a 2D array that contains y-direction
 *                   noise to be added to the translation. If provided, the
 *                   noise values are added to `dy` on a per-element basis.
 * @param  bbox      Domain bounding box.
 *
 * @return           A new 2D array that is the result of translating the input
 *                   `array` by the specified `dx` and `dy` values.
 *
 * **Example**
 * @include ex_translate.cpp
 *
 * **Result**
 * @image html ex_translate.png
 */
Array translate(const Array &array,
                float        dx,
                float        dy,
                bool         periodic = false,
                const Array *p_noise_x = nullptr,
                const Array *p_noise_y = nullptr,
                glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Apply a warping effect to the array.
 *
 * This function applies a warping effect to the input array by translating its
 * elements according to the specified x and y translation arrays. The warp
 * effect distorts the array based on the displacement values provided by `p_dx`
 * and `p_dy`.
 *
 * @param array Input array to be warped.
 * @param p_dx  Pointer to the array containing x-axis translation values.
 * @param p_dy  Pointer to the array containing y-axis translation values.
 *
 * **Example**
 * @include ex_warp.cpp
 *
 * **Result**
 * @image html ex_warp.png
 */
void warp(Array &array, const Array *p_dx, const Array *p_dy);

/**
 * @brief Apply a warping effect following the downward local gradient direction
 * (deflate/inflate effect).
 *
 * This function applies a warping effect to the input array based on the
 * downward local gradient direction, creating a deflate or inflate effect. The
 * warp amount, pre-filtering radius, and displacement direction can be
 * customized.
 *
 * @param array   Input array to be warped.
 * @param angle   The angle to determine the gradient direction.
 * @param amount  Amount of displacement (default is 0.02f).
 * @param ir      Pre-filtering radius to smooth the input data (default is 4).
 * @param reverse Reverse the displacement direction if set to true (default is
 *                false).
 *
 * **Example**
 * @include ex_warp_directional.cpp
 *
 * **Result**
 * @image html ex_warp_directional.png
 */
void warp_directional(Array &array,
                      float  angle,
                      float  amount = 0.02f,
                      int    ir = 4,
                      bool   reverse = false);

/**
 * @brief Apply a warping effect following the downward local gradient direction
 * (deflate/inflate effect) with a mask.
 *
 * This overloaded function applies a warping effect to the input array using a
 * specified mask. The mask controls the regions where the warp effect is
 * applied, with values expected in the range [0, 1]. The warp is based on the
 * downward local gradient direction.
 *
 * @param array   Input array to be warped.
 * @param angle   The angle to determine the gradient direction.
 * @param p_mask  Pointer to the mask array that filters the effect, expected in
 *                [0, 1].
 * @param amount  Amount of displacement (default is 0.02f).
 * @param ir      Pre-filtering radius to smooth the input data (default is 4).
 * @param reverse Reverse the displacement direction if set to true (default is
 *                false).
 *
 * **Example**
 * @include ex_warp_directional.cpp
 *
 * **Result**
 * @image html ex_warp_directional.png
 */
void warp_directional(Array       &array,
                      float        angle,
                      const Array *p_mask,
                      float        amount = 0.02f,
                      int          ir = 4,
                      bool         reverse = false); ///< @overload

/**
 * @brief Apply a warping effect following the downward local gradient direction
 * (deflate/inflate effect).
 *
 * This function applies a warping effect to the input array based on the
 * downward local gradient direction, simulating a deflate or inflate effect.
 * The effect can be customized by adjusting the displacement amount,
 * pre-filtering radius, and whether the displacement direction is reversed.
 *
 * @param array   Input array to be warped.
 * @param amount  Amount of displacement (default is 0.02f).
 * @param ir      Pre-filtering radius to smooth the input data (default is 4).
 * @param reverse Reverse the displacement direction if set to true (default is
 *                false).
 *
 * **Example**
 * @include ex_warp_downslope.cpp
 *
 * **Result**
 * @image html ex_warp_downslope.png
 */
void warp_downslope(Array &array,
                    float  amount = 0.02f,
                    int    ir = 4,
                    bool   reverse = false);

/**
 * @brief Apply a warping effect following the downward local gradient direction
 * (deflate/inflate effect) with a mask.
 *
 * This overloaded function applies a warping effect to the input array based on
 * the downward local gradient direction using a specified mask. The mask
 * controls where the warp effect is applied, with values expected in the range
 * [0, 1]. This function allows for additional customization of the warp effect.
 *
 * @param array   Input array to be warped.
 * @param p_mask  Pointer to the mask array that filters the effect, expected in
 *                [0, 1].
 * @param amount  Amount of displacement (default is 0.02f).
 * @param ir      Pre-filtering radius to smooth the input data (default is 4).
 * @param reverse Reverse the displacement direction if set to true (default is
 *                false).
 *
 * **Example**
 * @include ex_warp_downslope.cpp
 *
 * **Result**
 * @image html ex_warp_downslope.png
 */
void warp_downslope(Array       &array,
                    const Array *p_mask,
                    float        amount = 0.02f,
                    int          ir = 4,
                    bool         reverse = false); ///< @overload

/**
 * @brief Applies a zoom effect to a 2D array with an adjustable center.
 *
 * This function scales the input 2D array by a specified zoom factor,
 * effectively resizing the array's contents. The zoom operation is centered
 * around a specified point within the array, allowing for flexible zooming
 * behavior. The function supports both periodic boundary conditions, where the
 * array wraps around, and non-periodic conditions, where areas outside the
 * original array bounds are filled with zeros.
 *
 * @param  array       The input 2D array to be zoomed. This array remains
 *                     unmodified.
 * @param  zoom_factor The factor by which to zoom the array. A value greater
 *                     than 1 enlarges the contents, while a value between 0 and
 *                     1 reduces them.
 * @param  periodic    If set to `true`, the zoom is periodic, meaning that
 *                     elements moving out of the array bounds due to the zoom
 *                     reappear on the opposite side. If `false`, areas outside
 *                     the original array bounds are filled with zeros. The
 *                     default is `false`.
 * @param  center      The center of the zoom operation, specified as a
 *                     `glm::vec2`
 * with coordinates in the range [0, 1], where {0.5f, 0.5f} represents the
 * center of the array. The default center is {0.5f, 0.5f}.
 * @param  p_noise_x   Optional pointer to a 2D array that contains x-direction
 *                     noise to be added during the zoom operation. If provided,
 *                     the noise values are applied to the x-coordinates of the
 *                     zoomed array on a per-element basis.
 * @param  p_noise_y   Optional pointer to a 2D array that contains y-direction
 *                     noise to be added during the zoom operation. If provided,
 *                     the noise values are applied to the y-coordinates of the
 *                     zoomed array on a per-element basis.
 * @param  bbox        Domain bounding box.
 *
 * @return             A new 2D array that is the result of applying the zoom
 *                     effect to the input `array` by the specified
 *                     `zoom_factor` and centered at the specified
 * `center`.
 *
 * **Example**
 * @include ex_zoom.cpp
 *
 * **Result**
 * @image html ex_zoom.png
 */
Array zoom(const Array &array,
           float        zoom_factor,
           bool         periodic = false,
           glm::vec2    center = {0.5f, 0.5f},
           const Array *p_noise_x = nullptr,
           const Array *p_noise_y = nullptr,
           glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Performs particle-based advection on a scalar field.
 *
 * This function advects particles through a velocity field (or gradient field)
 * derived from the input array @p z. Each particle carries a value from the @p
 * advected_field and evolves over a number of iterations according to advection
 * parameters such as inertia and persistence. Optionally, advection and masking
 * fields can be applied to control which areas are affected.
 *
 * @param  z                 Input scalar field defining the base flow or
 *                           gradient field used for particle advection.
 * @param  advected_field    The scalar field whose values are transported along
 *                           particle trajectories.
 * @param  iterations        Number of advection iterations to perform per
 *                           particle.
 * @param  nparticles        Number of particles to simulate during advection.
 * @param  seed              Random seed used for initializing particle
 *                           positions and stochastic behavior.
 * @param  reverse           If true, reverses the advection direction (useful
 *                           for backtracing or inverse advection).
 * @param  post_filter       If true, applies a smoothing or filtering operation
 *                           after advection to reduce artifacts.
 * @param  advection_length  Maximum advection step length per iteration
 *                           (controls particle displacement scale).
 * @param  value_persistence Controls how much of a particle’s advected value is
 *                           preserved between iterations (1.0 means full
 *                           preservation, lower values introduce decay).
 * @param  inertia           Controls how much a particle’s previous velocity
 *                           affects its next step (simulates momentum).
 * @param  p_advection_mask  Optional mask array controlling where advection
 *                           occurs (nullptr to disable).
 * @param  p_mask            Optional general mask restricting where output
 *                           values are written (nullptr to disable).
 *
 * @return                   An Array representing the advected version of @p
 *                           advected_field after applying particle advection.
 *
 * **Example**
 * @include ex_advection_particle.cpp
 *
 * **Result**
 * @image html ex_advection_particle.png
 */
Array advection_particle(const Array  &z,
                         const Array  &advected_field,
                         int           iterations,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse = false,
                         bool          post_filter = true,
                         float         post_filter_sigma = 0.125f,
                         float         advection_length = 0.1f,
                         float         value_persistence = 0.99f,
                         float         inertia = 0.f,
                         const Array  *p_advection_mask = nullptr,
                         const Array  *p_mask = nullptr);

Array advection_particle(const Array  &z,
                         const Array  &advected_field,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse = false,
                         bool          post_filter = true,
                         float         post_filter_sigma = 0.125f,
                         float         advection_length = 0.1f,
                         float         value_persistence = 0.99f,
                         float         inertia = 0.f,
                         const Array  *p_advection_mask = nullptr,
                         const Array  *p_mask = nullptr);

Array advection_particle(const Array  &dx,
                         const Array  &dy,
                         const Array  &advected_field,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse = false,
                         bool          post_filter = true,
                         float         post_filter_sigma = 0.125f,
                         float         advection_length = 0.1f,
                         float         value_persistence = 0.99f,
                         float         inertia = 0.f,
                         const Array  *p_advection_mask = nullptr,
                         const Array  *p_mask = nullptr);

/**
 * @brief Performs 2D field advection based on the gradient of a heightmap using
 * a warp-based technic (simplified approach).
 *
 * This function computes the X and Y gradients of the given heightmap `z` and
 * uses them as flow directions to advect the input `advected_field`. The result
 * is a new field where values have been displaced according to the gradient
 * direction and the specified advection parameters.
 *
 * @param  z                 The heightmap (or scalar field) used to compute
 *                           advection directions.
 * @param  advected_field    The field to be advected along the gradient of `z`.
 * @param  advection_length  The step length used for advection; larger values
 *                           produce stronger displacement.
 * @param  value_persistence The persistence factor applied to the advected
 *                           values, typically in the range [0, 1]; controls how
 *                           much of the original value is preserved.
 * @param  p_mask            Optional pointer to a mask array; if provided,
 *                           advection is applied only where the mask is
 *                           nonzero.
 *
 * @return                   A new Array representing the advected version of
 *                           `advected_field`.
 *
 * **Example**
 * @include ex_advection_warp.cpp
 *
 * **Result**
 * @image html ex_advection_warp.png
 */
Array advection_warp(const Array &z,
                     const Array &advected_field,
                     float        advection_length = 0.1f,
                     float        value_persistence = 0.9f,
                     const Array *p_mask = nullptr);

Array advection_warp(const Array &z,
                     const Array &advected_field,
                     const Array &dx,
                     const Array &dy,
                     float        advection_length = 0.1f,
                     float        value_persistence = 0.9f,
                     const Array *p_mask = nullptr);

/*! @brief See hmap::rotate */
void rotate(Array &array, float angle, bool zoom_in = true);

/*! @brief See hmap::warp */
void warp(Array &array, const Array *p_dx, const Array *p_dy);

} // namespace hmap::gpu
