/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file rock_distribution.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstdint>
#include <string>

#include "highmap/rocks/rock.hpp"

namespace hmap
{

/**
 * @struct RockDistribution
 * @brief Parameters defining the size, power-law scaling, and embedding traits
 * of rocks.
 */
struct RockDistribution
{
  uint32_t    class_id = 0; ///< Rock material / mesh class identifier.
  std::string name =
      ""; ///< Optional descriptive name (e.g., "Granite Boulder").
  float radius_min = 1e-3f; ///< Minimum rock radius.
  float radius_max = 5e-2f; ///< Maximum rock radius.
  float power_law_alpha =
      2.0f; ///< Power-law / Pareto size exponent (typically in [1.5, 3.0]).
  float embed_ratio = 0.2f; ///< Ratio of rock embedded into the ground [0, 1].

  RockDistribution() = default;

  RockDistribution(uint32_t           class_id,
                   float              radius_min = 1e-3f,
                   float              radius_max = 5e-2f,
                   float              power_law_alpha = 2.0f,
                   float              embed_ratio = 0.2f,
                   const std::string &name = "")
      : class_id(class_id),
        name(name),
        radius_min(radius_min),
        radius_max(radius_max),
        power_law_alpha(power_law_alpha),
        embed_ratio(embed_ratio)
  {
  }
};

} // namespace hmap
