/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file rock_procedural.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstddef>

#include "highmap/array.hpp"
#include "highmap/rocks/rock_field.hpp"
#include "highmap/rocks/rock_seeding.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

// ============================================================================
//  Procedural Rock Placement Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Computes a probability density map in [0, 1] for rock placement based
 * on elevation, talus/slope, river drainage paths, an optional secondary
 * density map, and an optional exclusion mask.
 *
 * A subtle bulkify filter can be applied to the heightmap prior to computing
 * talus and river beds, eliminating micro-depressions and guiding flow lines.
 *
 * @param  heightmap            Input terrain heightmap.
 * @param  max_slope            Talus value above which density smoothly falls
 *                              to zero.
 * @param  bulkify_amp          Amplitude of bulkify pulse filter applied to the
 *                              terrain (0.f disables).
 * @param  weight_elev          Weight for inverted elevation criterion.
 * @param  weight_slope         Weight for talus/slope criterion.
 * @param  weight_rivers        Weight for river accumulation criterion.
 * @param  river_talus_ref      Reference talus for river flow partitioning.
 * @param  river_clipping_ratio Flow accumulation clipping ratio.
 * @param  secondary_density    Optional secondary density modulation map.
 * @param  weight_secondary     Weight for optional secondary density map.
 * @param  exclusion_mask       Optional exclusion mask (density set to 0 where
 *                              exclusion > exclusion_threshold).
 * @param  exclusion_threshold  Threshold above which exclusion mask excludes
 *                              rocks.
 * @param  bbox                 Domain bounding box for tiled bulkify synthesis.
 * @return                      Array Computed density array in [0, 1].
 */
Array build_rock_density(const Array &heightmap,
                         float        max_slope = 0.05f,
                         float        bulkify_amp = 0.01f,
                         float        weight_elev = 1.0f,
                         float        weight_slope = 1.0f,
                         float        weight_rivers = 1.0f,
                         float        river_talus_ref = 0.1f,
                         float        river_clipping_ratio = 1000.0f,
                         const Array &secondary_density = {},
                         float        weight_secondary = 1.0f,
                         const Array &exclusion_mask = {},
                         float        exclusion_threshold = 0.f,
                         glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Places rocks procedurally on a terrain, without physical simulation.
 *
 * Positions are sampled from a density map (defaults to a flat density equal to
 * 1 if none is provided). Rock radii are drawn from a bounded power-law
 * (Pareto) distribution between radius_min and a local maximum radius that
 * grows toward radius_max as the elevation decreases: high areas only hold
 * small rocks, whereas low areas hold both small and large rocks (the size
 * distribution broadens when lowering altitude). Steep slopes further reduce
 * the maximum rock size.
 *
 * @param  rock_count              Target number of rocks.
 * @param  elevation               Terrain elevation heightmap.
 * @param  min_elev                Elevation mapped to 0 in the normalized
 *                                 elevation.
 * @param  max_elev                Elevation mapped to 1 in the normalized
 *                                 elevation.
 * @param  max_slope               Talus threshold above which rock size is
 *                                 attenuated.
 * @param  density                 Optional user density map (if empty, defaults
 *                                 to uniform flat density of 1).
 * @param  exclusion               Optional exclusion map (same shape as
 *                                 elevation).
 * @param  size_slope_weight       In [0, 1], reduces the maximum rock size on
 *                                 steep slopes (0 = no effect).
 * @param  size_elevation_gain     Low-elevation influence on the maximum rock
 *                                 size (0 = none).
 * @param  size_elevation_exponent Power-law exponent of the maximum size /
 *                                 (1 - elevation) relation.
 * @param  options                 Seeding options (bbox, seed, exclusion
 *                                 threshold, rock distribution).
 * @return                         RockField Placed rocks.
 */
RockField place_rocks_procedural(size_t       rock_count,
                                 const Array &elevation,
                                 float        min_elev = 0.f,
                                 float        max_elev = 0.9f,
                                 float        max_slope = 1.f,
                                 const Array &density = {},
                                 const Array &exclusion = {},
                                 float        size_slope_weight = 0.f,
                                 float        size_elevation_gain = 1.f,
                                 float        size_elevation_exponent = 1.f,
                                 const RockSeedingOptions &options = {});

} // namespace hmap

namespace hmap::va
{

/**
 * @brief Places rocks procedurally across a tiled VirtualArray domain using
 * local density/elevation and bounding-box-aware tile merging.
 *
 * @param  rock_count              Target number of rocks.
 * @param  elevation               Terrain elevation VirtualArray.
 * @param  min_elev                Elevation mapped to 0 in normalized
 *                                 elevation.
 * @param  max_elev                Elevation mapped to 1 in normalized
 *                                 elevation.
 * @param  max_slope               Talus threshold above which rock size is
 *                                 attenuated.
 * @param  density                 Optional user density VirtualArray (if empty,
 *                                 defaults to uniform flat density of 1).
 * @param  exclusion               Optional exclusion VirtualArray (same shape
 *                                 as elevation).
 * @param  size_slope_weight       In [0, 1], reduces maximum rock size on steep
 *                                 slopes.
 * @param  size_elevation_gain     Low-elevation influence on maximum rock size.
 * @param  size_elevation_exponent Power-law exponent of size / (1 - elevation)
 *                                 relation.
 * @param  options                 Seeding options.
 * @param  cm                      Virtual array compute mode.
 * @return                         RockField Placed rocks.
 */
RockField place_rocks_procedural(size_t              rock_count,
                                 const VirtualArray &elevation,
                                 float               min_elev = 0.f,
                                 float               max_elev = 0.9f,
                                 float               max_slope = 1.f,
                                 const VirtualArray &density = {},
                                 const VirtualArray &exclusion = {},
                                 float               size_slope_weight = 0.f,
                                 float               size_elevation_gain = 1.f,
                                 float size_elevation_exponent = 1.f,
                                 const RockSeedingOptions &options = {},
                                 const ComputeMode        &cm = {});

} // namespace hmap::va
