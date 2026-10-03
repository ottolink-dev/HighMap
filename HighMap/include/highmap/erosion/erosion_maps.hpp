/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Computes erosion and deposition delta maps from terrain states before
 * and after erosion.
 *
 * @param z_before       Input array (before erosion).
 * @param z_after        Input array (after erosion).
 * @param erosion_map    Erosion map.
 * @param deposition_map Deposition map.
 * @param tolerance      Tolerance for erosion / deposition definition.
 *
 * **Example**
 * @include ex_erosions_maps.cpp
 *
 * **Result**
 * @image html ex_erosions_maps0.png
 * @image html ex_erosions_maps1.png
 */
void erosion_maps(Array &z_before,
                  Array &z_after,
                  Array &erosion_map,
                  Array &deposition_map,
                  float  tolerance = 0.f);

/**
 * @brief Generates a modified bedrock heightmap from an input elevation array.
 *
 * This function adjusts elevations based on overall height and local slope:
 * - Reduces values relative to the elevation range using `elevation_strength`.
 * - Reduces values in steep areas using `slope_strength` and `slope_limit`.
 *
 * @param  z                  Input elevation array.
 * @param  elevation_strength Strength of elevation-based adjustment.
 * @param  slope_strength     Strength of slope-based adjustment.
 * @param  slope_limit        Slope threshold for slope-based adjustment.
 * @param  zmin               Minimum elevation to consider (computed from z if
 *                            zmin > zmax).
 * @param  zmax               Maximum elevation to consider (computed from z if
 *                            zmin > zmax).
 * @return                    Array Modified bedrock heightmap.
 *
 *  **Example**
 * @include ex_hydraulic_particle.cpp
 *
 * **Result**
 * @image html ex_hydraulic_particle.png
 */
Array generate_bedrock(const Array &z,
                       float        elevation_strength,
                       float        slope_strength,
                       float        slope_limit,
                       float        zmin = 0.f,
                       float        zmax = -1.f);

} // namespace hmap
