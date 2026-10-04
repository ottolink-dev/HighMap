/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "highmap/geometry/inverse_sampler_2d.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/rocks/rock_seeding.hpp"
#include "highmap/scatter/scatter_seeding.hpp"

namespace hmap
{

// ============================================================================
//  Rock Seeding Functions (Alphabetically Sorted)
// ============================================================================

RockField seed_rock_field(size_t                    rock_count,
                          const Array              &density,
                          const Array              &exclusion,
                          const RockSeedingOptions &options)
{
  if (rock_count == 0 || !validate_non_empty(density)) return RockField();

  ScatterSeedingOptions scatter_opts;
  scatter_opts.bbox = options.bbox;
  scatter_opts.seed = options.seed;
  scatter_opts.exclusion_threshold = options.exclusion_threshold;
  scatter_opts.default_radius = options.distribution.radius_min;

  ScatterField field = seed_scatter_clusters(1,
                                             rock_count,
                                             density,
                                             exclusion,
                                             0.05f,
                                             8,
                                             scatter_opts);

  RockField rock_field(std::move(field));
  rock_field.apply_power_law_distribution(options.distribution, options.seed);
  rock_field.prune_collisions();

  return rock_field;
}

} // namespace hmap
