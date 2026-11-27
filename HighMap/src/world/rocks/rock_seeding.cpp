/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/scatter/scatter_field.hpp"
#include "highmap/scatter/scatter_seeding.hpp"
#include "highmap/virtual_array/virtual_array.hpp"
#include "highmap/world/rocks/rock_distribution.hpp"
#include "highmap/world/rocks/rock_field.hpp"
#include "highmap/world/rocks/rock_seeding.hpp"

namespace hmap
{

// ============================================================================
//  Rock Seeding Functions (Alphabetically Sorted)
// ============================================================================

RockField seed_rock_field(size_t                    rock_count,
                          const Array              &density,
                          const Array              &exclusion,
                          const RockSeedingOptions &options,
                          const UvTile             &region,
                          ICommunicator            &comm)
{
  if (rock_count == 0 || !validate_non_empty(density)) return RockField();

  float width = options.bbox.max.x - options.bbox.min.x;
  float height = options.bbox.max.y - options.bbox.min.y;
  if (width <= 1e-7f || height <= 1e-7f) return RockField();

  BBox2D local_bbox(options.bbox.min.x + region.uv_rect.u_min * width,
                    options.bbox.min.x + region.uv_rect.u_max * width,
                    options.bbox.min.y + region.uv_rect.v_min * height,
                    options.bbox.min.y + region.uv_rect.v_max * height);

  // --- Modulate density by exclusion map for collective quota calculation

  Array effective_density = density;
  bool  has_exclusion = !exclusion.vector.empty() || exclusion.shape.x > 0;
  if (has_exclusion)
  {
    if (!validate_non_empty(exclusion) ||
        !validate_same_shape(density, exclusion))
      return RockField();

    for (int j = 0; j < effective_density.shape.y; ++j)
      for (int i = 0; i < effective_density.shape.x; ++i)
      {
        if (exclusion(i, j) > options.exclusion_threshold)
          effective_density(i, j) = 0.0f;
      }
  }

  // --- Collective density reduction & local rock quota

  float local_sum = effective_density.sum();
  float total_sum = comm.all_reduce(local_sum, ReduceOp::Sum);
  if (total_sum <= 1e-7f) return RockField();

  size_t local_rock_count = rock_count;
  if (comm.size() > 1)
  {
    local_rock_count = static_cast<size_t>(
        std::round(static_cast<float>(rock_count) * (local_sum / total_sum)));
  }
  if (local_rock_count == 0) return RockField();

  uint32_t local_seed = (comm.size() > 1)
                            ? options.seed + static_cast<uint32_t>(comm.rank())
                            : options.seed;

  ScatterSeedingOptions scatter_opts;
  scatter_opts.bbox = local_bbox;
  scatter_opts.seed = local_seed;
  scatter_opts.exclusion_threshold = options.exclusion_threshold;
  scatter_opts.default_radius = options.distribution.radius_min;

  ScatterField field = seed_scatter_clusters(1,
                                             local_rock_count,
                                             density,
                                             exclusion,
                                             0.05f,
                                             8,
                                             scatter_opts);

  RockField rock_field(std::move(field));
  rock_field.apply_power_law_distribution(options.distribution, local_seed);
  rock_field.prune_collisions();

  return rock_field;
}

} // namespace hmap
