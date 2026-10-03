/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include "highmap/array.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

/**
 * @brief Fill the depressions of the heightmap using the Planchon-Darboux
 * algorithm.
 *
 * Fill heightmap depressions to ensure that every cell can be connected to the
 * boundaries following a downward slope @cite Planchon2002.
 *
 * @param z              Input array.
 * @param iterations     Number of iterations.
 * @param epsilon        Epsilon slope parameter.
 * @param outflow_left   Allow outflow on left boundary (i = 0).
 * @param outflow_right  Allow outflow on right boundary (i = nx - 1).
 * @param outflow_bottom Allow outflow on bottom boundary (j = 0).
 * @param outflow_top    Allow outflow on top boundary (j = ny - 1).
 *
 * **Example**
 * @include ex_depression_filling.cpp
 *
 * **Result**
 * @image html ex_depression_filling.png
 */
void depression_filling(Array &z,
                        int    iterations = 1000,
                        float  epsilon = 1e-4f,
                        bool   outflow_left = true,
                        bool   outflow_right = true,
                        bool   outflow_bottom = true,
                        bool   outflow_top = true);

/**
 * @brief Fill depressions in a heightmap using the Priority-Flood algorithm.
 *
 * @param z                 Input heightmap array.
 * @param apply_post_filter Apply Laplacian smoothing to deposition.
 * @param outflow_left      Allow outflow on left boundary (i = 0).
 * @param outflow_right     Allow outflow on right boundary (i = nx - 1).
 * @param outflow_bottom    Allow outflow on bottom boundary (j = 0).
 * @param outflow_top       Allow outflow on top boundary (j = ny - 1).
 */
void depression_filling_priority_flood(Array &z,
                                       bool   apply_post_filter = false,
                                       bool   outflow_left = true,
                                       bool   outflow_right = true,
                                       bool   outflow_bottom = true,
                                       bool   outflow_top = true);

} // namespace hmap

namespace hmap::va
{

/**
 * @brief Fill depressions in a VirtualArray using the Priority-Flood algorithm.
 *
 * @param  z                 Input heightmap.
 * @param  apply_post_filter Apply Laplacian smoothing to deposition.
 * @param  p_fill_map        Optional output fill map (z_after - z_before).
 * @param  cm                Compute mode configuration.
 * @return                   Filled heightmap.
 */
VirtualArray depression_filling_priority_flood(
    const VirtualArray &z,
    bool                apply_post_filter = false,
    VirtualArray       *p_fill_map = nullptr,
    const ComputeMode  &cm = {});

} // namespace hmap::va
