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

namespace hmap::va
{

RockField seed_rock_field(size_t                    rock_count,
                          const VirtualArray       &density,
                          const VirtualArray       &exclusion,
                          const RockSeedingOptions &options,
                          const ComputeMode        &cm)
{
  // --- Validate inputs

  if (rock_count == 0 || density.empty())
  {
    return RockField();
  }

  bool has_exclusion = (exclusion.storage != nullptr) &&
                       (exclusion.shape.x > 0 && exclusion.shape.y > 0);

  if (has_exclusion)
  {
    if (exclusion.shape != density.shape ||
        exclusion.tile_shape != density.tile_shape)
    {
      hmap::log::warn(
          "Exclusion VirtualArray shape mismatch with density VirtualArray");
      return RockField();
    }
  }

  glm::ivec2 num_tiles = density.get_max_tiles();
  int        total_tiles = num_tiles.x * num_tiles.y;

  if (total_tiles == 0) return RockField();

  // --- Compute Density Sums per Tile and Global Inner Sum

  const VirtualArray *p_exclusion_va = has_exclusion ? &exclusion : nullptr;

  std::vector<float> tile_full_sums;
  float              total_inner_sum = 0.f;

  sum_per_tile(density,
               tile_full_sums,
               &total_inner_sum,
               p_exclusion_va,
               options.exclusion_threshold,
               cm);

  if (total_inner_sum <= 1e-7f) return RockField();

  // --- Generate Sub-RockFields per Tile

  std::vector<ScatterField> tile_fields(total_tiles);
  std::vector<glm::vec4>    tile_bboxs(total_tiles);

  hmap::for_each_tile(
      {&density, p_exclusion_va},
      {},
      [&](const std::vector<const Array *> &p_arrays_in,
          std::vector<Array *> &,
          const UvTile &region)
      {
        const auto [p_density_tile, p_exclusion_tile] = unpack<2>(p_arrays_in);

        int tile_index = region.key.ty * num_tiles.x + region.key.tx;

        tile_bboxs[tile_index] = region.uv_rect.to_vec4();

        float  tile_sum = tile_full_sums[tile_index];
        size_t tile_rock_count = static_cast<size_t>(std::round(
            static_cast<float>(rock_count) * (tile_sum / total_inner_sum)));

        if (tile_rock_count == 0)
        {
          return;
        }

        RockSeedingOptions tile_options = options;
        tile_options.bbox = region.uv_rect.to_vec4();
        tile_options.seed = options.seed +
                            static_cast<uint32_t>(tile_index) * 1999u + 17u;

        RockField tile_rock_field = seed_rock_field(
            tile_rock_count,
            *p_density_tile,
            p_exclusion_tile ? *p_exclusion_tile : Array(),
            tile_options);

        tile_fields[tile_index] = std::move(tile_rock_field);
      },
      cm);

  // --- Merge RockFields with Bounding Box Overlap Resolution

  ScatterField merged = merge_scatter_fields(tile_fields, tile_bboxs, true);

  return RockField(std::move(merged));
}

} // namespace hmap::va
