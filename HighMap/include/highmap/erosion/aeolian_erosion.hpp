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
 * @param z                   Input/output total elevation heightmap (modified
 *                            in-place). On output, @p z equals bedrock +
 *                            sediment.
 * @param steps               Number of simulation iterations.
 * @param seed                Random seed for stochastic sampling.
 * @param wind                Base wind direction vector (in grid units).
 * @param p_bedrock           Optional input bedrock elevation array. If null,
 *                            bedrock is initialized to zero.
 * @param p_vegetation        Optional vegetation density mask in [0, 1].
 * @param p_abrasion_mask     Optional rock erodibility/susceptibility mask in
 *                            [0, 1].
 * @param p_sediment_map      Optional output array receiving final sand
 *                            thickness.
 * @param p_bedrock_map       Optional output array receiving final bedrock
 *                            elevation.
 * @param cell_size           Horizontal grid spacing in world units. If <= 0,
 *                            defaults to 1.f / z.shape.x.
 * @param matter_to_move      Volume/height of sand grain packet displaced per
 *                            hop.
 * @param samples_per_step    Number of stochastic saltation samples evaluated
 *                            per step. If <= 0, defaults to total cell count
 *                            (nx *
 *                            ny).
 * @param sand_repose_angle   Angle of repose for loose sediment (in degrees,
 *                            typically ~33°).
 * @param shadow_angle_min    Minimum wind shadow angle threshold (in degrees,
 *                            typically ~5°).
 * @param shadow_angle_max    Maximum wind shadow angle threshold (in degrees,
 *                            typically ~15°).
 * @param max_shadow_dist     Maximum raymarch distance for wind shadow
 *                            detection (in grid cells).
 * @param abrasion_rate       Abrasion wear rate on exposed bedrock (0 to
 *                            disable).
 * @param repose_rate         Rate of sediment redistribution during repose
 *                            stabilization.
 * @param sync_interval       Number of simulation steps between communicator
 *                            tile syncs.
 * @param region              Tile region geometry metadata.
 * @param comm                Tiled communicator instance.
 * @param shadow_map_interval If > 0, wind shadows are precomputed on a grid
 *                            every @p shadow_map_interval steps and looked up
 *                            by the walkers, instead of ray-marching per walker
 *                            hop. Much faster, at the cost of a slightly
 *                            approximate (stale and interpolated) shadow. 0
 *                            keeps the exact per-walker ray-marching.
 * @param repose_interval     Number of simulation steps between angle-of-repose
 *                            stabilization passes.
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
                     ICommunicator    &comm = SingleArrayComm::instance(),
                     int               shadow_map_interval = 1,
                     int               repose_interval = 2);

/**
 * @brief Multiscale (coarse-to-fine) variant of hmap::aeolian_erosion.
 *
 * The simulation is first run on downsampled copies of the terrain, where each
 * step is much cheaper (fewer cells and shorter shadow rays in cells), to
 * develop the large-scale dune pattern. At each level, only the *change* in
 * sediment (and bedrock) is upsampled and added to the full-resolution layers
 * so that fine details of the input are preserved. The last steps are run at
 * full resolution to restore small-scale features.
 *
 * @param z                   Input/output total elevation (bedrock + sediment).
 * @param steps               Total number of steps at full resolution (steps
 *                            are distributed across the levels).
 * @param seed                Random seed.
 * @param wind                Base wind vector (reference grid of 256).
 * @param p_bedrock           Optional bedrock elevation.
 * @param p_vegetation        Optional vegetation density in [0, 1].
 * @param p_abrasion_mask     Optional abrasion susceptibility in [0, 1].
 * @param p_sediment_map      Optional output sediment thickness.
 * @param p_bedrock_map       Optional output bedrock elevation.
 * @param cell_size           Full-resolution cell size (<= 0 for 1 /
 *                            z.shape.x).
 * @param matter_to_move      Sand packet displaced per hop.
 * @param sand_repose_angle   Angle of repose (degrees).
 * @param shadow_angle_min    Minimum wind shadow angle (degrees).
 * @param shadow_angle_max    Maximum wind shadow angle (degrees).
 * @param max_shadow_dist     Maximum shadow raymarch distance (reference grid).
 * @param abrasion_rate       Bedrock abrasion rate (0 to disable).
 * @param repose_rate         Repose stabilization rate.
 * @param levels              Number of resolution levels (1 = single-scale).
 *                            Each level halves the resolution.
 * @param coarse_fraction     Fraction in [0, 1) of the steps spent on coarse
 *                            levels (shared equally); the remainder is run at
 *                            full resolution.
 * @param min_level_size      Smallest allowed size of a coarse level (levels
 *                            that would be smaller are dropped).
 * @param shadow_map_interval See hmap::aeolian_erosion (0 = exact shadows).
 * @param repose_interval     See hmap::aeolian_erosion.
 */
void aeolian_erosion_multiscale(Array        &z,
                                int           steps = 300,
                                std::uint32_t seed = 1337u,
                                glm::vec2     wind = {3.f, 0.f},
                                const Array  *p_bedrock = nullptr,
                                const Array  *p_vegetation = nullptr,
                                const Array  *p_abrasion_mask = nullptr,
                                Array        *p_sediment_map = nullptr,
                                Array        *p_bedrock_map = nullptr,
                                float         cell_size = -1.f,
                                float         matter_to_move = 0.001f,
                                float         sand_repose_angle = 33.f,
                                float         shadow_angle_min = 5.f,
                                float         shadow_angle_max = 15.f,
                                float         max_shadow_dist = 12.f,
                                float         abrasion_rate = 0.f,
                                float         repose_rate = 0.35f,
                                int           levels = 3,
                                float         coarse_fraction = 0.7f,
                                int           min_level_size = 64,
                                int           shadow_map_interval = 1,
                                int           repose_interval = 2);

} // namespace hmap
