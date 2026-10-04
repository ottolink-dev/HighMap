/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file rock_seeding.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/rocks/rock_distribution.hpp"
#include "highmap/rocks/rock_field.hpp"

namespace hmap
{

/**
 * @struct RockSeedingOptions
 * @brief Configuration parameters for rock field seeding algorithms.
 */
struct RockSeedingOptions
{
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}; ///< Domain bounding box.
  uint32_t  seed = 0;                    ///< Random number generator seed.
  float     exclusion_threshold = 0.5f;  ///< Threshold for exclusion map masking.
  RockDistribution distribution =
      {}; ///< Rock size & power-law distribution traits.
};

// ============================================================================
//  Rock Seeding Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Seeds a general rock/boulder field using inverse density transform
 * sampling and applies power-law radius distribution and collision resolution.
 *
 * @param  rock_count Target number of rocks.
 * @param  density    Spatial probability density map array.
 * @param  exclusion  Optional exclusion map array.
 * @param  options    Seeding options.
 * @return            RockField Sampled rock field.
 */
RockField seed_rock_field(size_t                    rock_count,
                          const Array              &density,
                          const Array              &exclusion = {},
                          const RockSeedingOptions &options = {});

} // namespace hmap
