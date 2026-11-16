/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file species.hpp
 * @copyright Copyright (c) 2025 Otto Link.
 */
#pragma once

#include <cstdint>
#include <string>

#include "highmap/flora/tree.hpp"

namespace hmap
{

/**
 * @struct Species
 * @brief Encapsulates biological and geometric traits of a tree species.
 */
struct Species
{
  uint32_t    id = 0;    ///< Species identifier.
  std::string name = ""; ///< Optional human-readable name.
  float       radius =
      HMAP_DEFAULT_TREE_RADIUS; ///< Nominal canopy / collision radius.
  float weight = 1.0f;          ///< Relative abundance weight for seeding.
  float radius_min = 0.5f *
                     HMAP_DEFAULT_TREE_RADIUS; ///< Minimum viable crown radius.
  float radius_max = 1.5f * HMAP_DEFAULT_TREE_RADIUS; ///< Maximum crown radius.
  float competition_factor = 0.4f; ///< Default intra/inter-species alpha.

  Species() = default;

  Species(uint32_t           id,
          float              radius = HMAP_DEFAULT_TREE_RADIUS,
          float              weight = 1.0f,
          float              radius_min = -1.0f,
          float              radius_max = -1.0f,
          float              competition_factor = 0.4f,
          const std::string &name = "")
      : id(id),
        name(name),
        radius(radius),
        weight(weight),
        radius_min((radius_min >= 0.0f) ? radius_min : 0.5f * radius),
        radius_max((radius_max >= 0.0f) ? radius_max : 1.5f * radius),
        competition_factor(competition_factor)
  {
  }
};

} // namespace hmap
