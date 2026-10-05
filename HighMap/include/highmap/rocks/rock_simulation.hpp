/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file rock_simulation.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/rocks/rock_distribution.hpp"

namespace hmap
{

class RockField;

/**
 * @struct RockSimulationOptions
 * @brief Configuration parameters for GPU-based rock and debris trajectory
 * simulation.
 */
struct RockSimulationOptions
{
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};    ///< Domain bounding box {xmin, xmax, ymin, ymax}.
  float     time_step = 0.005f;             ///< Simulation time step dt in seconds.
  int       max_steps = 1500;               ///< Maximum total simulation steps.
  float     gravity = 9.81f;                ///< Gravitational acceleration.
  float     soil_friction = 0.55f;          ///< Coulomb friction coefficient mu.
  float     rolling_resistance = 0.05f;     ///< Rolling resistance coefficient.
  float     inter_rock_restitution = 0.20f; ///< Restitution between colliding rocks.
  float min_velocity = 0.01f;               ///< Resting threshold speed.
  bool  respawn_out_of_bounds = false;      ///< Whether to respawn out-of-bounds rocks at initial position.
  float spawn_fraction = 0.0f;              ///< Fraction of max_steps over which rocks
                               ///< progressively spawn (0.0 = all at start).
  uint32_t seed = 0; ///< Random seed.
};

// ============================================================================
//  Rock Simulation Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Emits rocks from a spatial emission density map and simulates their
 * trajectories downhill across terrain on GPU.
 *
 * @param  rock_count       Target number of rocks to emit and simulate.
 * @param  emission_density Spatial probability density map where rocks
 *                          originate.
 * @param  elevation        Terrain elevation heightmap.
 * @param  p_friction_map   Optional soil friction map (same dimensions as
 *                          elevation).
 * @param  dist             Rock distribution specifying radii and classes.
 * @param  options          Simulation options.
 * @return                  RockField       Field containing the resting rocks.
 */
RockField simulate_rock_emission(size_t       rock_count,
                                 const Array &emission_density,
                                 const Array &elevation,
                                 const Array *p_friction_map = nullptr,
                                 const RockDistribution      &dist = {},
                                 const RockSimulationOptions &options = {});

/**
 * @brief Simulates physical movement (sliding, rolling, bouncing, and
 * inter-rock collision avoidance) for an existing RockField across terrain on
 * GPU.
 *
 * @param  rocks          Input rock field with initial positions.
 * @param  elevation      Terrain elevation heightmap.
 * @param  p_friction_map Optional soil friction map.
 * @param  options        Simulation options.
 * @return                RockField     Updated rock field with resting
 *                        positions.
 */
RockField simulate_rock_trajectories(const RockField &rocks,
                                     const Array     &elevation,
                                     const Array     *p_friction_map = nullptr,
                                     const RockSimulationOptions &options = {});

} // namespace hmap
