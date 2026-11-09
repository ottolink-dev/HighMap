/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <cstdint>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/virtual_array/tile_region.hpp"

namespace hmap
{

class ICommunicator;

/**
 * @brief Simulates wind-driven aeolian terrain erosion, sand transport, and
 * dune dynamics on GPU.
 *
 * This algorithm is adapted from the research by Axel Paris, Adrien Peytavie,
 * Eric Guérin, Oscar Argudo, and Eric Galin (Pacific Graphics 2019 / Computer
 * Graphics Forum), "Desertscapes Simulation".
 *
 * It simulates:
 * - Wind-directed sand saltation and reptation across the landscape.
 * - Lee-side aerodynamic wind shadows triggering dune deposition.
 * - Vegetation interaction, anchoring sand to form vegetated mounds (nabkhas).
 * - Bedrock abrasion in high-energy, sediment-sparse regions, carving yardangs.
 * - Angle-of-repose avalanching maintaining realistic slipfaces.
 *
 * @param z                 Input/output total elevation heightmap (modified
 *                          in-place). On output, @p z equals bedrock + sediment.
 * @param steps             Number of simulation iterations.
 * @param seed              Random seed for stochastic sampling.
 * @param wind              Base wind direction vector (in grid units).
 * @param p_bedrock         Optional input bedrock elevation array. If null,
 *                          bedrock is initialized to zero.
 * @param p_vegetation      Optional vegetation density mask in [0, 1].
 * @param p_abrasion_mask   Optional rock erodibility/susceptibility mask in [0, 1].
 * @param p_sediment_map    Optional output array receiving final sand thickness.
 * @param p_bedrock_map     Optional output array receiving final bedrock elevation.
 * @param cell_size         Horizontal grid spacing in world units. If <= 0,
 *                          defaults to 1.f / z.shape.x.
 * @param matter_to_move    Volume/height of sand grain packet displaced per hop.
 * @param samples_per_step  Number of stochastic saltation samples evaluated per step.
 *                          If <= 0, defaults to total cell count (nx * ny).
 * @param sand_repose_angle Angle of repose for loose sediment (in degrees, typically ~33°).
 * @param shadow_angle_min  Minimum wind shadow angle threshold (in degrees, typically ~5°).
 * @param shadow_angle_max  Maximum wind shadow angle threshold (in degrees, typically ~15°).
 * @param max_shadow_dist   Maximum raymarch distance for wind shadow detection (in grid cells).
 * @param abrasion_rate     Abrasion wear rate on exposed bedrock (0 to disable).
 * @param repose_rate       Rate of sediment redistribution during repose stabilization.
 * @param sync_interval     Number of simulation steps between communicator tile syncs.
 * @param region            Tile region geometry metadata.
 * @param comm              Tiled communicator instance.
 *
 * **Example**
 * @include ex_aeolian_erosion.cpp
 *
 * **Result**
 * @image html ex_aeolian_erosion.png
 */
void aeolian_erosion(Array            &z,
                     int               steps = 300,
                     std::uint32_t     seed = 1337u,
                     glm::vec2         wind = {3.f, 0.f},
                     const Array      *p_bedrock = nullptr,
                     const Array      *p_vegetation = nullptr,
                     const Array      *p_abrasion_mask = nullptr,
                     Array            *p_sediment_map = nullptr,
                     Array            *p_bedrock_map = nullptr,
                     float             cell_size = -1.f,
                     float             matter_to_move = 0.001f,
                     int               samples_per_step = 0,
                     float             sand_repose_angle = 33.f,
                     float             shadow_angle_min = 5.f,
                     float             shadow_angle_max = 15.f,
                     float             max_shadow_dist = 12.f,
                     float             abrasion_rate = 0.f,
                     float             repose_rate = 0.35f,
                     int               sync_interval = 0,
                     const TileRegion &region = TileRegion::dummy(),
                     ICommunicator    &comm = SingleArrayComm::instance());

} // namespace hmap
