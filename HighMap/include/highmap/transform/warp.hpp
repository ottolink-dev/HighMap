/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file warp.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for array warping functions.
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{

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

/// @overload
void warp_directional(Array       &array,
                      float        angle,
                      const Array *p_mask,
                      float        amount = 0.02f,
                      int          ir = 4,
                      bool         reverse = false);

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

/// @overload
void warp_downslope(Array       &array,
                    const Array *p_mask,
                    float        amount = 0.02f,
                    int          ir = 4,
                    bool         reverse = false);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::warp */
void warp(Array &array, const Array *p_dx, const Array *p_dy);

} // namespace hmap::gpu
