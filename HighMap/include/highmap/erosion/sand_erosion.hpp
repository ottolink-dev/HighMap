/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <cmath>
#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/erosion/erosion_parameters.hpp"

namespace hmap
{

/**
 * @brief Parameters for the sand deposition and cellular dune simulation model
 * (Werner 1995 extension with elevation-aware hop length and boundary control).
 */
struct SandDuneParams
{
  float wind_angle = 0.f;            // wind direction angle in radians (0 = along +x)
  float wind_speed = 1.0f;           // base wind speed factor
  float crest_speedup = 0.5f;        // speedup factor on ridges / higher elevations
  float hop_length = 5.0f;           // base transport hop distance in grid units
  float elevation_hop_factor = 2.0f; // extra hop distance per unit elevation drop downwind
  float jitter_angle = 0.2f;         // lateral stochastic wind jitter angle in radians
  float p_sand = 0.6f;               // deposition probability on sandy substrate
  float p_bare = 0.4f;               // deposition probability on bare bedrock
  float shadow_talus = 0.002f;       // aerodynamic shadow talus (height difference between 2 cells)
  float talus = 0.005f;              // angle of repose talus limit (height difference between 2 cells)
  float slab_height = 0.005f;        // erosion slab thickness per particle
  float collapse_rate = 0.5f;        // avalanche relaxation rate factor
  int   max_hops = 10;               // maximum consecutive saltation hops if not deposited
  bool  periodic_boundary = true;    // toroidal boundary wrapping (true) or respawn (false)
  float sand_inflow_rate = 0.0f;     // sand replenishment rate per pass at upwind border
  std::uint32_t seed = 42;           // random seed
};

namespace gpu
{

/**
 * @brief Simulates sand deposition and dune morphology (Werner cellular model)
 * on GPU.
 *
 * @param z                Heightmap array modified in-place.
 * @param nparticles       Total number of particle transport events.
 * @param params           Solver parameters.
 * @param p_bedrock        Optional bedrock array bounding maximum erosion depth.
 * @param p_mask           Optional blend / active area mask in [0, 1].
 * @param p_deposition_map Optional output array receiving cumulative deposition.
 * @param p_erosion_map    Optional output array receiving cumulative erosion.
 * @param iterations       Number of simulation passes (splitting particles across passes).
 */
void sand_dune(Array                &z,
               int                   nparticles = 50000,
               const SandDuneParams &params = SandDuneParams(),
               const Array          *p_bedrock = nullptr,
               const Array          *p_mask = nullptr,
               Array                *p_deposition_map = nullptr,
               Array                *p_erosion_map = nullptr,
               int                   iterations = 10);

} // namespace gpu

} // namespace hmap
