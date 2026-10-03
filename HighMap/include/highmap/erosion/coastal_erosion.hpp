/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Simulates terrain diffusion due to coastal erosion.
 *
 * This function applies an iterative coastal erosion diffusion process on a
 * terrain elevation array (`z`), taking into account the presence and depth of
 * water. The erosion model smooths the terrain near shorelines while
 * maintaining constant water-surface height.
 *
 * At each iteration:
 *  - A local water mask is computed using `water_mask(water_depth, z,
 * additional_depth)`.
 *  - The terrain elevations are smoothed (diffused) using a masked Laplacian
 * filter, applied only near the water boundary.
 *  - The water depth is adjusted to preserve the total water surface height
 * (i.e., `z + water_depth` remains constant).
 *
 * This results in a realistic simulation of coastal erosion processes where the
 * terrain near the waterline is progressively smoothed and redistributed.
 *
 * @param z                Reference to the terrain elevation array (modified in
 *                         place).
 * @param water_depth      Reference to the array representing water depth
 *                         values (modified to preserve water surface level).
 * @param additional_depth Additional virtual water depth used to estimate the
 *                         influence region of the shoreline during mask
 *                         computation.
 * @param iterations       Number of erosion–diffusion iterations to apply.
 * @param p_mask           Optional intensity mask in [0, 1].
 * @param p_water_mask     Optional output pointer. If non-null, it receives the
 *                         last computed water mask used during the final
 *                         iteration.
 *
 * **Example**
 * @include ex_coastal_erosion_diffusion.cpp
 *
 * **Result**
 * @image html ex_coastal_erosion_diffusion.png
 */
void coastal_erosion_diffusion(Array       &z,
                               Array       &water_depth,
                               float        additional_depth,
                               int          iterations = 10,
                               const Array *p_mask = nullptr,
                               Array       *p_water_mask = nullptr);

/**
 * @brief Applies a coastal erosion profile to a terrain elevation field.
 *
 * This function modifies the elevation array @p z by carving a coastal profile
 * at the interface between ground and water. It uses distance transforms both
 * from ground regions and from water regions to determine how far each point is
 * from the shoreline. According to this distance, the function applies:
 *
 * - A ground-side shore slope with an optional scarp transition.
 * - A water-side underwater slope ensuring continuity with the ground slope.
 * - An optional post-filtering step (Laplace smoothing) restricted to the
 * shoreline region.
 *
 * Water depth is adjusted after filtering to preserve the original water
 * surface height.
 *
 * Optionally, a shoreline mask can be returned describing the region influenced
 * by the coastal transformation.
 *
 * @param z                           Array of terrain elevations. Modified
 *                                    in-place.
 * @param water_depth                 Array of water depths. Modified to
 *                                    preserve water surface height after
 *                                    terrain changes.
 * @param shore_ground_extent         Horizontal extent (in grid units) over
 *                                    which the ground-side shore profile is
 *                                    applied.
 * @param shore_water_extent          Horizontal extent (in grid units) over
 *                                    which the underwater profile is applied.
 * @param slope_shore                 Ground-side slope magnitude of the coastal
 *                                    profile. Expressed in elevation units per
 *                                    domain width.
 * @param slope_shore_water           Water-side slope magnitude of the
 *                                    underwater profile. Expressed similarly to
 *                                    @p slope_shore.
 * @param scarp_extent_ratio          Ratio defining the relative extent of the
 *                                    scarp region. A value in [0,1]:
 *                                    - 0 → no scarp, only slope
 *                                    - 1 → all scarp
 * @param apply_post_filter           If true, applies Laplacian smoothing to @p
 *                                    z restricted to shoreline areas.
 * @param post_filter_iterations      Iterations for post-smoothing filter.
 * @param solid_shore_mask            Whether shore mask is solid.
 * @param scarp_mask_transition_ratio Smooth transition ratio for scarp mask.
 * @param p_noise                     Optional noise array.
 * @param p_shore_mask                Optional output pointer. If non-null,
 *                                    receives the shoreline mask (values in
 *                                    [0,1]) indicating where the coastal
 *                                    transformation was applied.
 * @param p_scarp_mask                Optional output pointer for scarp mask.
 *
 * **Example**
 * @include ex_coastal_erosion_profile.cpp
 *
 * **Result**
 * @image html ex_coastal_erosion_profile.png
 */
void coastal_erosion_profile(Array &z,
                             Array &water_depth,
                             float  shore_ground_extent, // pixels
                             float  shore_water_extent,
                             float  slope_shore = 0.5f,
                             float  slope_shore_water = 0.5f,
                             float  scarp_extent_ratio = 0.5f, // in [0,
                                                               // 1]
                             bool         apply_post_filter = true,
                             int          post_filter_iterations = 3,
                             bool         solid_shore_mask = true,
                             float        scarp_mask_transition_ratio = 0.2f,
                             const Array *p_noise = nullptr,
                             Array       *p_shore_mask = nullptr,
                             Array       *p_scarp_mask = nullptr);

void coastal_erosion_profile(Array       &z,
                             const Array *p_mask,
                             Array       &water_depth,
                             float        shore_ground_extent, // pixels
                             float        shore_water_extent,
                             float        slope_shore = 0.5f,
                             float        slope_shore_water = 0.5f,
                             float        scarp_extent_ratio = 0.5f, // in [0,
                                                                     // 1]
                             bool         apply_post_filter = true,
                             int          post_filter_iterations = 3,
                             bool         solid_shore_mask = true,
                             float        scarp_mask_transition_ratio = 0.2f,
                             const Array *p_noise = nullptr,
                             Array       *p_shore_mask = nullptr,
                             Array       *p_scarp_mask = nullptr);

} // namespace hmap
