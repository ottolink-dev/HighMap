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
  TBM_HEIGHT,
  TBM_MAX,
  TBM_MAX_SMOOTH,
  TBM_MIN,
  TBM_MIN_SMOOTH,
  TBM_MULTISCALE,
  TBM_POISSON,
  TBM_SEAMLESS,
  TBM_SLOPE,
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

/**
 * @brief Height-aware terrain blending for Array inputs.
 */
void blend_terrains_height(Array            &lower,
                           const CoordFrame &frame_lower,
                           const Array      &upper,
                           const CoordFrame &frame_upper,
                           float             feather_width = 0.2f,
                           float             height_offset = 0.f,
                           float             height_range = 0.2f);

/**
 * @brief Height-aware terrain blending for VirtualArray inputs.
 */
void blend_terrains_height(VirtualArray       &lower,
                           const CoordFrame   &frame_lower,
                           const VirtualArray &upper,
                           const CoordFrame   &frame_upper,
                           float               feather_width = 0.2f,
                           float               height_offset = 0.f,
                           float               height_range = 0.2f,
                           const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Multiscale (frequency-band) terrain blending for Array inputs.
 */
void blend_terrains_multiscale(Array            &lower,
                               const CoordFrame &frame_lower,
                               const Array      &upper,
                               const CoordFrame &frame_upper,
                               float             feather_width = 0.25f,
                               int               n_levels = 3,
                               int               base_radius = 16,
                               float             detail_gain = 1.f);

/**
 * @brief Multiscale (frequency-band) terrain blending for VirtualArray inputs
 * using a coarse global proxy for low-frequency bands and tile-local detail
 * blending with halo synchronization.
 */
void blend_terrains_multiscale(VirtualArray       &lower,
                               const CoordFrame   &frame_lower,
                               const VirtualArray &upper,
                               const CoordFrame   &frame_upper,
                               float               feather_width = 0.25f,
                               int                 n_levels = 3,
                               int                 base_radius = 16,
                               float               detail_gain = 1.f,
                               int                 coarse_size = 128,
                               const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Poisson / gradient-domain terrain blending for Array inputs.
 */
void blend_terrains_poisson(Array            &lower,
                            const CoordFrame &frame_lower,
                            const Array      &upper,
                            const CoordFrame &frame_upper,
                            int               iterations = 400,
                            float             feather_width = 0.15f,
                            bool              mixed_gradients = true);

/**
 * @brief Poisson / gradient-domain terrain blending for VirtualArray inputs
 * using a two-level coarse-proxy + tile-local Red-Black SOR solver with halo
 * synchronization.
 */
void blend_terrains_poisson(VirtualArray       &lower,
                            const CoordFrame   &frame_lower,
                            const VirtualArray &upper,
                            const CoordFrame   &frame_upper,
                            int                 iterations = 400,
                            float               feather_width = 0.15f,
                            bool                mixed_gradients = true,
                            int                 coarse_size = 128,
                            int                 sync_interval = 25,
                            const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Seamless boundary blending (harmonic MVC membrane) for Array inputs.
 */
void blend_terrains_seamless(Array            &lower,
                             const CoordFrame &frame_lower,
                             const Array      &upper,
                             const CoordFrame &frame_upper,
                             float             feather_width = 0.2f,
                             int               boundary_samples = 64);

/**
 * @brief Seamless boundary blending (harmonic MVC membrane) for VirtualArray
 * inputs.
 */
void blend_terrains_seamless(VirtualArray       &lower,
                             const CoordFrame   &frame_lower,
                             const VirtualArray &upper,
                             const CoordFrame   &frame_upper,
                             float               feather_width = 0.2f,
                             int                 boundary_samples = 64,
                             const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Slope-aware terrain blending for Array inputs.
 */
void blend_terrains_slope(Array            &lower,
                          const CoordFrame &frame_lower,
                          const Array      &upper,
                          const CoordFrame &frame_upper,
                          float             feather_width = 0.2f,
                          float             slope_scale = 0.1f,
                          int               ir_slope = 2);

/**
 * @brief Slope-aware terrain blending for VirtualArray inputs.
 */
void blend_terrains_slope(VirtualArray       &lower,
                          const CoordFrame   &frame_lower,
                          const VirtualArray &upper,
                          const CoordFrame   &frame_upper,
                          float               feather_width = 0.2f,
                          float               slope_scale = 0.1f,
                          int                 ir_slope = 2,
                          const ComputeMode  &cm = ComputeMode{});

} // namespace hmap
