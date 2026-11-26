/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file transform.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for array transformation functions including rotation,
 * flipping, symmetry, translation, transposition, and zooming.
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include <cstdint>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @enum SymmetryType
 * @brief Describes the symmetry operation to apply to an array.
 */
// clang-format off
enum SymmetryType : int
{
	SYMMETRY_LEFT_TO_RIGHT, ///< Mirror left half onto right half.
	SYMMETRY_RIGHT_TO_LEFT, ///< Mirror right half onto left half.
	SYMMETRY_TOP_TO_BOTTOM, ///< Mirror top half onto bottom half.
	SYMMETRY_BOTTOM_TO_TOP, ///< Mirror bottom half onto top half.
	SYMMETRY_X,             ///< Symmetric average across vertical axis.
	SYMMETRY_Y,             ///< Symmetric average across horizontal axis.
	SYMMETRY_XY,            ///< 4-quadrant symmetric average across both X and Y.
	SYMMETRY_ROT180         ///< 180-degree rotational symmetry around center.
};
// clang-format on

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
void radial_displacement_to_xy(const Array  &dr,
                               Array        &dx,
                               Array        &dy,
                               float         smoothing = 1.f,
                               glm::vec2     center = {0.5f, 0.5f},
                               const UvRect &uv_rect = UvRect::unit());

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
 * @brief Scale the array in UV space with wrapping boundary.
 *
 * @param array    Input array to be scaled.
 * @param uv_scale Scaling factors for U and V dimensions.
 */
void scale_uv(Array &array, glm::vec2 uv_scale);

/**
 * @brief Enforces geometric symmetry on an array.
 *
 * This function reflects or averages portions of the input array across the
 * specified symmetry axis or plane according to @p symmetry_type.
 *
 * @param  array          Input array to symmetrize.
 * @param  symmetry_type  Type of symmetry operation to apply.
 * @param  flatten_center If true, flattens the elevation along the symmetry
 *                        axis or center point towards minimum elevation.
 * @param  flatten_radius Radius (in normalized [0, 1] units) around the
 *                        symmetry axis/center for the flattening transition.
 * @return                Symmetrized array.
 */
Array symmetrize(
    const Array &array,
    SymmetryType symmetry_type = SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
    bool         flatten_center = false,
    float        flatten_radius = 0.05f);

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
Array translate(const Array  &array,
                float         dx,
                float         dy,
                bool          periodic = false,
                const Array  *p_noise_x = nullptr,
                const Array  *p_noise_y = nullptr,
                const UvRect &uv_rect = UvRect::unit());

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
 *                     with coordinates in the range [0, 1], where {0.5f, 0.5f}
 *                     represents the center of the array. The default center is
 *                     {0.5f, 0.5f}.
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
 *                     `center`.
 *
 * **Example**
 * @include ex_zoom.cpp
 *
 * **Result**
 * @image html ex_zoom.png
 */
Array zoom(const Array  &array,
           float         zoom_factor,
           bool          periodic = false,
           glm::vec2     center = {0.5f, 0.5f},
           const Array  *p_noise_x = nullptr,
           const Array  *p_noise_y = nullptr,
           const UvRect &uv_rect = UvRect::unit());

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::rotate */
void rotate(Array &array, float angle, bool zoom_in = true);

} // namespace hmap::gpu
