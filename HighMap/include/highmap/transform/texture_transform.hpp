/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file texture_transform.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for texture transformation functions.
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include <glm/vec2.hpp>

namespace hmap
{

class Texture;

/**
 * @brief Flip the texture horizontally (left/right).
 *
 * @param texture Input texture to be flipped.
 */
void flip_lr(Texture &texture);

/**
 * @brief Flip the texture vertically (up/down).
 *
 * @param texture Input texture to be flipped.
 */
void flip_ud(Texture &texture);

/**
 * @brief Rotate the texture by 180 degrees.
 *
 * @param texture Input texture to be rotated.
 */
void rot180(Texture &texture);

/**
 * @brief Rotate the texture by 270 degrees counterclockwise.
 *
 * @param texture Input texture to be rotated.
 */
void rot270(Texture &texture);

/**
 * @brief Rotate the texture by 90 degrees counterclockwise.
 *
 * @param texture Input texture to be rotated.
 */
void rot90(Texture &texture);

/**
 * @brief Scale the texture in UV space with wrapping boundary.
 *
 * @param texture  Input texture to be scaled.
 * @param uv_scale Scaling factors for U and V dimensions.
 */
void scale_uv(Texture &texture, glm::vec2 uv_scale);

/**
 * @brief Return the transposed texture.
 *
 * @param  texture Input texture to be transposed.
 * @return         Texture The transposed texture.
 */
Texture transpose(const Texture &texture);

} // namespace hmap
