/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#pragma once
#include "highmap/array.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/virtual_array.hpp"

namespace hmap
{

/**
 * @brief Blending modes for combining terrains with CoordFrame transforms.
 */
enum class TerrainBlendMode : int
{
  TBM_ADD,
  TBM_AVERAGE,
  TBM_DISTANCE,
  TBM_MAX,
  TBM_MAX_SMOOTH,
  TBM_MIN,
  TBM_MIN_SMOOTH,
  TBM_SUBTRACT,
};

/**
 * @brief Blend an upper terrain onto a lower terrain in-place, taking their
 * CoordFrame transforms into account.
 *
 * @param lower       Lower (base) terrain array, modified in-place.
 * @param frame_lower Coordinate frame of the lower terrain.
 * @param upper       Upper (overlay) terrain array.
 * @param frame_upper Coordinate frame of the upper terrain.
 * @param mode        Blending mode.
 * @param param       Mode-specific parameter (smoothness k for MIN_SMOOTH /
 *                    MAX_SMOOTH, normalized feathering width in (0, 0.5] for
 *                    DISTANCE).
 */
void blend_terrains(Array            &lower,
                    const CoordFrame &frame_lower,
                    const Array      &upper,
                    const CoordFrame &frame_upper,
                    TerrainBlendMode  mode,
                    float             param = 0.2f);

/**
 * @brief Blend an upper VirtualArray terrain onto a lower VirtualArray terrain
 * in-place, taking their CoordFrame transforms into account.
 *
 * @param lower       Lower (base) VirtualArray, modified in-place.
 * @param frame_lower Coordinate frame of the lower terrain.
 * @param upper       Upper (overlay) VirtualArray.
 * @param frame_upper Coordinate frame of the upper terrain.
 * @param mode        Blending mode.
 * @param param       Mode-specific parameter (smoothness k for MIN_SMOOTH /
 *                    MAX_SMOOTH, normalized feathering width in (0, 0.5] for
 *                    DISTANCE).
 * @param cm          Compute mode for tile execution.
 */
void blend_terrains(VirtualArray       &lower,
                    const CoordFrame   &frame_lower,
                    const VirtualArray &upper,
                    const CoordFrame   &frame_upper,
                    TerrainBlendMode    mode,
                    float               param = 0.2f,
                    const ComputeMode  &cm = ComputeMode{});

} // namespace hmap
