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
 * @brief Alignment modes for fitting an upper terrain to a lower terrain.
 */
enum class TerrainAlignMode : int
{
  TAM_OFFSET,
  TAM_PLANE,
  TAM_ROBUST_OFFSET,
  TAM_ROBUST_PLANE,
};

/**
 * @brief Elevation offset and tilt parameters for aligning an upper terrain.
 */
struct TerrainAlignment
{
  float offset = 0.f;
  float slope_x = 0.f;
  float slope_y = 0.f;
};

/**
 * @brief Compute alignment parameters to fit @p upper to @p lower over their
 * common valid area.
 */
TerrainAlignment compute_terrain_alignment(const Array      &lower,
                                           const CoordFrame &frame_lower,
                                           const Array      &upper,
                                           const CoordFrame &frame_upper,
                                           TerrainAlignMode  mode,
                                           float             huber_delta = 0.1f,
                                           int               max_iterations = 10);

/**
 * @brief Compute alignment parameters for VirtualArray inputs using a coarse
 * Array representation.
 */
TerrainAlignment compute_terrain_alignment(const VirtualArray &lower,
                                           const CoordFrame   &frame_lower,
                                           const VirtualArray &upper,
                                           const CoordFrame   &frame_upper,
                                           TerrainAlignMode    mode,
                                           float               huber_delta = 0.1f,
                                           int                 max_iterations = 10,
                                           int                 coarse_size = 128);

/**
 * @brief Apply precomputed alignment parameters to an Array in-place.
 */
void apply_terrain_alignment(Array                  &upper,
                             const TerrainAlignment &alignment);

/**
 * @brief Apply precomputed alignment parameters to a VirtualArray in-place.
 */
void apply_terrain_alignment(VirtualArray           &upper,
                             const TerrainAlignment &alignment,
                             const ComputeMode      &cm = ComputeMode{});

/**
 * @brief Align @p upper to @p lower in-place over their common valid area.
 */
void align_terrains(const Array      &lower,
                    const CoordFrame &frame_lower,
                    Array            &upper,
                    const CoordFrame &frame_upper,
                    TerrainAlignMode  mode,
                    float             huber_delta = 0.1f,
                    int               max_iterations = 10);

/**
 * @brief Align @p upper VirtualArray to @p lower VirtualArray in-place.
 */
void align_terrains(const VirtualArray &lower,
                    const CoordFrame   &frame_lower,
                    VirtualArray       &upper,
                    const CoordFrame   &frame_upper,
                    TerrainAlignMode    mode,
                    float               huber_delta = 0.1f,
                    int                 max_iterations = 10,
                    int                 coarse_size = 128,
                    const ComputeMode  &cm = ComputeMode{});

} // namespace hmap
