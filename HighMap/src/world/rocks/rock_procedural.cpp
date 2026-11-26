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

Array build_rock_density(const Array  &heightmap,
                         float         max_slope,
                         float         bulkify_amp,
                         float         weight_elev,
                         float         weight_slope,
                         float         weight_rivers,
                         float         river_talus_ref,
                         float         river_clipping_ratio,
                         const Array  &secondary_density,
                         float         weight_secondary,
                         const Array  &exclusion_mask,
                         float         exclusion_threshold,
                         const BBox2D &bbox)
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
                     UvRect::from_vec4(bbox.to_vec4()));
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
                                 const RockSeedingOptions &options,
                                 const UvTile             &region,
                                 ICommunicator            &comm)
{
  if (rock_count == 0 || !validate_non_empty(elevation)) return RockField();

  float width = options.bbox.max.x - options.bbox.min.x;
  float height = options.bbox.max.y - options.bbox.min.y;
  if (width <= 1e-7f || height <= 1e-7f) return RockField();

  BBox2D local_bbox(options.bbox.min.x + region.uv_rect.u_min * width,
                    options.bbox.min.x + region.uv_rect.u_max * width,
                    options.bbox.min.y + region.uv_rect.v_min * height,
                    options.bbox.min.y + region.uv_rect.v_max * height);

  // --- Density map: use user density if supplied, otherwise flat density equal
  // to 1

  Array target_density;
  bool  has_density = !density.vector.empty() || density.shape.x > 0;
  if (has_density)
  {
    if (!validate_non_empty(density) ||
        !validate_same_shape(elevation, density))
      return RockField();
    target_density = density;
  }
  else
  {
    target_density = Array(elevation.shape, 1.0f);
  }

  // --- Modulate density by exclusion map for collective quota calculation

  Array effective_density = target_density;
  bool  has_exclusion = !exclusion.vector.empty() || exclusion.shape.x > 0;
  if (has_exclusion)
  {
    if (!validate_non_empty(exclusion) ||
        !validate_same_shape(elevation, exclusion))
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

  // --- Normalized elevation and slope

  Array z_norm = threshold_smooth(elevation, min_elev, max_elev);

  Array slope = gradient_norm(elevation);
  slope = threshold_smooth(slope, 0.f, max_slope);

  // --- Sample positions

  ScatterSeedingOptions scatter_opts;
  scatter_opts.bbox = local_bbox;
  scatter_opts.seed = local_seed;
  scatter_opts.exclusion_threshold = options.exclusion_threshold;
  scatter_opts.default_radius = options.distribution.radius_min;

  ScatterField field = seed_scatter_clusters(1,
                                             local_rock_count,
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

  std::mt19937                          gen(local_seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  const UvRect local_uv_rect = UvRect::from_vec4(local_bbox.to_vec4());

  for (auto &rock : rock_field)
  {
    const float e = z_norm.get_value_bilinear(rock.position.x,
                                              rock.position.y,
                                              local_uv_rect);

    const float t = std::clamp(size_elevation_gain *
                                   std::pow(1.f - e, size_elevation_exponent),
                               0.f,
                               1.f);

    // steep slopes reduce the maximum rock size
    const float s = slope.get_value_bilinear(rock.position.x,
                                             rock.position.y,
                                             local_uv_rect);

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
                                                   local_uv_rect);
  }

  return rock_field;
}

} // namespace hmap
