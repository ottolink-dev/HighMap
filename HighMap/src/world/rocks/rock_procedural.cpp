/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/range.hpp"
#include "highmap/scatter/scatter_field.hpp"
#include "highmap/scatter/scatter_seeding.hpp"
#include "highmap/selector.hpp"
#include "highmap/virtual_array/virtual_array.hpp"
#include "highmap/world/rocks/rock_field.hpp"
#include "highmap/world/rocks/rock_procedural.hpp"
#include "highmap/world/rocks/rock_seeding.hpp"

namespace hmap
{

// ============================================================================
//  Procedural Rock Placement Functions (Alphabetically Sorted)
// ============================================================================

Array build_rock_density(const Array &heightmap,
                         float        max_slope,
                         float        bulkify_amp,
                         float        weight_elev,
                         float        weight_slope,
                         float        weight_rivers,
                         float        river_talus_ref,
                         float        river_clipping_ratio,
                         const Array &secondary_density,
                         float        weight_secondary,
                         const Array &exclusion_mask,
                         float        exclusion_threshold,
                         glm::vec4    bbox)
{
  if (!validate_non_empty(heightmap)) return Array();

  // --- Bulkify Filter (to eliminate micro-depressions & guide flow)

  Array z_work;
  if (bulkify_amp > 0.f)
  {
    z_work = bulkify(heightmap,
                     PrimitiveType::PRIM_CUBIC_PULSE,
                     bulkify_amp,
                     nullptr,
                     nullptr,
                     {0.5f, 0.5f},
                     bbox);
  }
  else
  {
    z_work = heightmap;
  }

  // --- Elevation Criterion (prefer valleys & low elevation)

  Array cz;
  if (weight_elev > 0.f)
  {
    cz = -z_work;
    remap(cz);
  }

  // --- Slope Criterion (penalize steep slopes)

  Array cg;
  if (weight_slope > 0.f)
  {
    Array slope = gradient_norm(z_work);
    cg = 1.0f - threshold_smooth(slope, 0.f, max_slope);
  }

  // --- River Bed / Drainage Criterion

  Array cr;
  if (weight_rivers > 0.f)
  {
    cr = select_rivers(z_work, river_talus_ref, river_clipping_ratio);
    remap(cr);
  }

  // --- Combine Density Layers

  std::vector<const Array *> layers;
  std::vector<float>         weights;

  if (weight_elev > 0.f)
  {
    layers.push_back(&cz);
    weights.push_back(weight_elev);
  }

  if (weight_slope > 0.f)
  {
    layers.push_back(&cg);
    weights.push_back(weight_slope);
  }

  if (weight_rivers > 0.f)
  {
    layers.push_back(&cr);
    weights.push_back(weight_rivers);
  }

  bool has_secondary = !secondary_density.vector.empty() ||
                       secondary_density.shape.x > 0;
  if (has_secondary && weight_secondary > 0.f)
  {
    if (!validate_non_empty(secondary_density) ||
        !validate_same_shape(heightmap, secondary_density))
      return Array();

    layers.push_back(&secondary_density);
    weights.push_back(weight_secondary);
  }

  Array density;
  if (layers.empty())
  {
    density = Array(heightmap.shape, 1.0f);
  }
  else
  {
    density = build_density_log(layers, weights);
  }

  // --- Exclusion Mask

  bool has_exclusion = !exclusion_mask.vector.empty() ||
                       exclusion_mask.shape.x > 0;
  if (has_exclusion)
  {
    if (!validate_non_empty(exclusion_mask) ||
        !validate_same_shape(heightmap, exclusion_mask))
      return Array();

    for (size_t i = 0; i < density.vector.size(); ++i)
    {
      if (exclusion_mask.vector[i] > exclusion_threshold)
      {
        density.vector[i] = 0.f;
      }
    }
  }

  return density;
}

RockField place_rocks_procedural(size_t       rock_count,
                                 const Array &elevation,
                                 float        min_elev,
                                 float        max_elev,
                                 float        max_slope,
                                 const Array &density,
                                 const Array &exclusion,
                                 float        size_slope_weight,
                                 float        size_elevation_gain,
                                 float        size_elevation_exponent,
                                 const RockSeedingOptions &options)
{
  if (rock_count == 0 || !validate_non_empty(elevation)) return RockField();

  // --- Density map: use user density if supplied, otherwise flat density equal
  // to 1

  Array target_density;
  bool  has_density = !density.vector.empty() || density.shape.x > 0;
  if (has_density)
  {
    target_density = density;
  }
  else
  {
    target_density = Array(elevation.shape, 1.0f);
  }

  // --- Normalized elevation and slope

  Array z_norm = threshold_smooth(elevation, min_elev, max_elev);

  Array slope = gradient_norm(elevation);
  slope = threshold_smooth(slope, 0.f, max_slope);

  // --- Sample positions

  ScatterSeedingOptions scatter_opts;
  scatter_opts.bbox = options.bbox;
  scatter_opts.seed = options.seed;
  scatter_opts.exclusion_threshold = 0.f;
  scatter_opts.default_radius = options.distribution.radius_min;

  ScatterField field = seed_scatter_clusters(1,
                                             rock_count,
                                             target_density,
                                             exclusion,
                                             0.05f,
                                             8,
                                             scatter_opts);

  RockField rock_field(std::move(field));

  // --- Elevation-driven radius: bounded power-law whose upper bound
  // --- grows when the elevation decreases

  const float radius_min = std::max(options.distribution.radius_min, 1e-9f);
  const float radius_max = std::max(radius_min,
                                    options.distribution.radius_max);
  const float alpha = std::max(options.distribution.power_law_alpha, 1e-3f);

  std::mt19937                          gen(options.seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  for (auto &rock : rock_field)
  {
    const float e = z_norm.get_value_bilinear(rock.position.x,
                                              rock.position.y,
                                              options.bbox);

    const float t = std::clamp(size_elevation_gain *
                                   std::pow(1.f - e, size_elevation_exponent),
                               0.f,
                               1.f);

    // steep slopes reduce the maximum rock size
    const float s = slope.get_value_bilinear(rock.position.x,
                                             rock.position.y,
                                             options.bbox);

    const float slope_factor = 1.f -
                               std::clamp(size_slope_weight, 0.f, 1.f) * s;

    const float local_radius_max = radius_min +
                                   t * slope_factor * (radius_max - radius_min);

    // inverse transform sampling of P(R >= r) ~ r^-alpha on
    // [radius_min, local_radius_max]
    const float a_min = std::pow(radius_min, -alpha);
    const float a_max = std::pow(local_radius_max, -alpha);
    const float radius = std::pow(a_min - dis(gen) * (a_min - a_max),
                                  -1.f / alpha);

    rock.radius = std::clamp(radius, radius_min, local_radius_max);
    rock.class_id = options.distribution.class_id;
    rock.position.z = elevation.get_value_bilinear(rock.position.x,
                                                   rock.position.y,
                                                   options.bbox);
  }

  return rock_field;
}

} // namespace hmap

namespace hmap::va
{

RockField place_rocks_procedural(size_t              rock_count,
                                 const VirtualArray &elevation,
                                 float               min_elev,
                                 float               max_elev,
                                 float               max_slope,
                                 const VirtualArray &density,
                                 const VirtualArray &exclusion,
                                 float               size_slope_weight,
                                 float               size_elevation_gain,
                                 float               size_elevation_exponent,
                                 const RockSeedingOptions &options,
                                 const ComputeMode        &cm)
{
  // --- Validate inputs

  if (rock_count == 0 || elevation.empty())
  {
    return RockField();
  }

  bool has_density = (density.storage != nullptr) &&
                     (density.shape.x > 0 && density.shape.y > 0);

  if (has_density)
  {
    if (density.shape != elevation.shape ||
        density.tile_shape != elevation.tile_shape)
    {
      hmap::log::warn(
          "Density VirtualArray shape mismatch with elevation VirtualArray");
      return RockField();
    }
  }

  bool has_exclusion = (exclusion.storage != nullptr) &&
                       (exclusion.shape.x > 0 && exclusion.shape.y > 0);

  if (has_exclusion)
  {
    if (exclusion.shape != elevation.shape ||
        exclusion.tile_shape != elevation.tile_shape)
    {
      hmap::log::warn(
          "Exclusion VirtualArray shape mismatch with elevation VirtualArray");
      return RockField();
    }
  }

  glm::ivec2 num_tiles = elevation.get_max_tiles();
  int        total_tiles = num_tiles.x * num_tiles.y;

  if (total_tiles == 0) return RockField();

  // --- Compute Density Sums per Tile and Global Inner Sum

  const VirtualArray *p_density_va = has_density ? &density : nullptr;
  const VirtualArray *p_exclusion_va = has_exclusion ? &exclusion : nullptr;

  std::vector<float> tile_full_sums;
  float              total_inner_sum = 0.f;

  if (has_density)
  {
    sum_per_tile(density,
                 tile_full_sums,
                 &total_inner_sum,
                 p_exclusion_va,
                 options.exclusion_threshold,
                 cm);
  }
  else
  {
    // Uniform density: count valid cells per tile
    tile_full_sums.assign(total_tiles, 0.f);

    hmap::for_each_tile(
        {&elevation, p_exclusion_va},
        {},
        [&](const std::vector<const Array *> &p_arrays_in,
            std::vector<Array *> &,
            const TileRegion &region)
        {
          const auto [p_elev_tile, p_excl_tile] = unpack<2>(p_arrays_in);

          int tile_index = region.key.ty * num_tiles.x + region.key.tx;

          float tile_full_sum = 0.f;
          float tile_inner_sum = 0.f;

          const glm::ivec4 &halo = region.halo;
          int               shape_x = region.shape.x;
          int               shape_y = region.shape.y;

          for (int j = 0; j < shape_y; ++j)
          {
            bool is_inner_y = (j >= halo.z && j < shape_y - halo.w);

            for (int i = 0; i < shape_x; ++i)
            {
              float val = 1.0f;

              if (p_excl_tile &&
                  (*p_excl_tile)(i, j) > options.exclusion_threshold)
              {
                val = 0.0f;
              }

              if (val > 0.0f)
              {
                tile_full_sum += val;

                if (is_inner_y && i >= halo.x && i < shape_x - halo.y)
                {
                  tile_inner_sum += val;
                }
              }
            }
          }

          tile_full_sums[tile_index] = tile_full_sum;
          total_inner_sum += tile_inner_sum;
        },
        cm);
  }

  if (total_inner_sum <= 1e-7f) return RockField();

  // --- Generate Sub-RockFields per Tile

  std::vector<ScatterField> tile_fields(total_tiles);
  std::vector<glm::vec4>    tile_bboxs(total_tiles);

  hmap::for_each_tile(
      {&elevation, p_density_va, p_exclusion_va},
      {},
      [&](const std::vector<const Array *> &p_arrays_in,
          std::vector<Array *> &,
          const TileRegion &region)
      {
        const auto [p_elev_tile, p_density_tile, p_excl_tile] = unpack<3>(
            p_arrays_in);

        int tile_index = region.key.ty * num_tiles.x + region.key.tx;

        tile_bboxs[tile_index] = region.bbox;

        float  tile_sum = tile_full_sums[tile_index];
        size_t tile_rock_count = static_cast<size_t>(std::round(
            static_cast<float>(rock_count) * (tile_sum / total_inner_sum)));

        if (tile_rock_count == 0)
        {
          return;
        }

        RockSeedingOptions tile_options = options;
        tile_options.bbox = region.bbox;
        tile_options.seed = options.seed +
                            static_cast<uint32_t>(tile_index) * 1999u + 17u;

        RockField tile_rock_field = hmap::place_rocks_procedural(
            tile_rock_count,
            *p_elev_tile,
            min_elev,
            max_elev,
            max_slope,
            p_density_tile ? *p_density_tile : Array(),
            p_excl_tile ? *p_excl_tile : Array(),
            size_slope_weight,
            size_elevation_gain,
            size_elevation_exponent,
            tile_options);

        tile_fields[tile_index] = std::move(tile_rock_field);
      },
      cm);

  // --- Merge RockFields with Bounding Box Overlap Resolution

  ScatterField merged = merge_scatter_fields(tile_fields, tile_bboxs, true);

  return RockField(std::move(merged));
}

} // namespace hmap::va
