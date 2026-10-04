/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>

#include "point_sampler/kmeans_clustering.hpp"
#include "point_sampler/metrics.hpp"
#include "point_sampler/utils.hpp"

#include "highmap/geometry/inverse_sampler_2d.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/scatter/scatter_seeding.hpp"

namespace hmap
{

// ============================================================================
//  Seeding Functions (Alphabetically Sorted)
// ============================================================================

ScatterField seed_scatter_clusters(size_t       class_count,
                                   size_t       item_count,
                                   const Array &density,
                                   const Array &exclusion,
                                   float        cluster_spread,
                                   size_t       points_per_cluster,
                                   const ScatterSeedingOptions &options)
{
  // validate inputs
  if (class_count == 0 || item_count == 0 || points_per_cluster == 0 ||
      !validate_non_empty(density))
    return ScatterField();

  float width = options.bbox.y - options.bbox.x;
  float height = options.bbox.w - options.bbox.z;
  if (width <= 1e-7f || height <= 1e-7f) return ScatterField();

  // --- Modulate Density by Exclusion Map

  Array effective_density = density;
  bool  has_exclusion = !exclusion.vector.empty() || exclusion.shape.x > 0;

  if (has_exclusion)
  {
    if (!validate_non_empty(exclusion) ||
        !validate_same_shape(density, exclusion))
      return ScatterField();

    for (int j = 0; j < density.shape.y; ++j)
      for (int i = 0; i < density.shape.x; ++i)
      {
        if (exclusion(i, j) >= options.exclusion_threshold)
          effective_density(i, j) = 0.0f;
      }
  }

  InverseSampler2D sampler(effective_density, options.seed, options.bbox);
  if (sampler.get_total_weight() <= 1e-7f) return ScatterField();

  // --- Determine Class Quotas

  std::vector<size_t> class_item_counts(class_count, 0);
  if (options.class_weights.size() == class_count)
  {
    float total_w = 0.0f;
    for (float w : options.class_weights)
      total_w += std::max(0.0f, w);

    if (total_w > 1e-6f)
    {
      size_t allocated = 0;
      for (size_t s = 0; s < class_count; ++s)
      {
        float norm_w = std::max(0.0f, options.class_weights[s]) / total_w;
        class_item_counts[s] = static_cast<size_t>(
            std::round(static_cast<float>(item_count) * norm_w));
        allocated += class_item_counts[s];
      }

      // adjust difference due to rounding
      if (allocated < item_count)
        class_item_counts[0] += (item_count - allocated);
      else if (allocated > item_count &&
               class_item_counts[0] >= (allocated - item_count))
        class_item_counts[0] -= (allocated - item_count);
    }
  }
  else
  {
    size_t base_count = item_count / class_count;
    size_t remainder = item_count % class_count;
    for (size_t s = 0; s < class_count; ++s)
      class_item_counts[s] = base_count + (s < remainder ? 1 : 0);
  }

  // --- Generate Clusters per Class

  std::vector<ScatterItem> candidate_items;
  candidate_items.reserve(item_count);

  for (size_t s = 0; s < class_count; ++s)
  {
    if (class_item_counts[s] == 0) continue;

    size_t target_s = class_item_counts[s];
    size_t cluster_count = std::max<size_t>(
        1,
        (target_s + points_per_cluster - 1) / points_per_cluster);

    auto parent_samples = sampler.sample(cluster_count);

    size_t count_for_s = 0;
    for (size_t c = 0; c < cluster_count && count_for_s < target_s; ++c)
    {
      float cx = parent_samples[0][c];
      float cy = parent_samples[1][c];

      size_t pts_to_sample = std::min(points_per_cluster,
                                      target_s - count_for_s);
      if (pts_to_sample == 0) break;

      glm::vec4 cluster_bbox = {
          cx - cluster_spread,
          cx + cluster_spread,
          cy - cluster_spread,
          cy + cluster_spread,
      };

      for (size_t k = 0; k < pts_to_sample; ++k)
      {
        glm::vec2 pt = sampler.sample(cluster_bbox);
        candidate_items.emplace_back(pt.x,
                                     pt.y,
                                     0.0f,
                                     static_cast<uint32_t>(s),
                                     options.default_radius);
        count_for_s++;
      }
    }
  }

  return ScatterField(std::move(candidate_items));
}

ScatterField seed_scatter_kmeans(size_t       class_count,
                                 size_t       item_count,
                                 const Array &density,
                                 const Array &exclusion,
                                 float        cluster_randomness,
                                 size_t       k_neighbors,
                                 const ScatterSeedingOptions &options)
{
  // validate inputs
  if (class_count == 0 || item_count == 0 || !validate_non_empty(density))
    return ScatterField();

  float width = options.bbox.y - options.bbox.x;
  float height = options.bbox.w - options.bbox.z;
  if (width <= 1e-7f || height <= 1e-7f) return ScatterField();

  // --- Modulate Density by Exclusion Map

  Array effective_density = density;
  bool  has_exclusion = !exclusion.vector.empty() || exclusion.shape.x > 0;

  if (has_exclusion)
  {
    if (!validate_non_empty(exclusion) ||
        !validate_same_shape(density, exclusion))
      return ScatterField();

    for (int j = 0; j < density.shape.y; ++j)
      for (int i = 0; i < density.shape.x; ++i)
      {
        if (exclusion(i, j) >= options.exclusion_threshold)
          effective_density(i, j) = 0.0f;
      }
  }

  InverseSampler2D sampler(effective_density, options.seed, options.bbox);
  if (sampler.get_total_weight() <= 1e-7f) return ScatterField();

  // --- Sample Coordinates Using 2D Inverse Sampling

  auto   samples = sampler.sample(item_count);
  size_t actual_count = samples[0].size();
  if (actual_count == 0) return ScatterField();

  // --- Cluster and Tag Classes

  std::vector<ScatterItem> candidate_items;
  candidate_items.reserve(actual_count);

  if (class_count == 1 || actual_count == 1)
  {
    for (size_t i = 0; i < actual_count; ++i)
    {
      candidate_items.emplace_back(samples[0][i],
                                   samples[1][i],
                                   0.0f,
                                   0u,
                                   options.default_radius);
    }
  }
  else
  {
    // construct point list for spatial neighbor analysis
    std::vector<ps::Point<float, 2>> points;
    points.reserve(actual_count);
    for (size_t i = 0; i < actual_count; ++i)
    {
      points.push_back(ps::Point<float, 2>({samples[0][i], samples[1][i]}));
    }

    size_t safe_k_neighbors = std::min(std::max<size_t>(1, k_neighbors),
                                       actual_count - 1);

    auto idx = ps::nearest_neighbors_indices(points, safe_k_neighbors);

    std::vector<float> dist_min;
    std::vector<float> dist_avg;
    dist_min.reserve(idx.size());
    dist_avg.reserve(idx.size());

    float min_dmin = 1e9f, max_dmin = 0.f;
    float min_davg = 1e9f, max_davg = 0.f;

    for (size_t k = 0; k < idx.size(); ++k)
    {
      float  dmin = 1e9f;
      float  davg = 0.f;
      size_t n_nbrs = idx[k].size();

      for (size_t r = 0; r < n_nbrs; ++r)
      {
        float dist = ps::distance_squared(points[k], points[idx[k][r]]);
        dmin = std::min(dist, dmin);
        davg += dist;
      }

      dist_min.push_back(dmin);
      dist_avg.push_back(davg);

      min_dmin = std::min(min_dmin, dmin);
      max_dmin = std::max(max_dmin, dmin);
      min_davg = std::min(min_davg, davg);
      max_davg = std::max(max_davg, davg);
    }

    // add scalable random component in [0, 1] to features if requested
    float rand_amount = std::clamp(cluster_randomness, 0.0f, 1.0f);

    if (rand_amount > 0.0f)
    {
      std::mt19937                          rng(options.seed + 1);
      std::uniform_real_distribution<float> dist_uni(-0.5f, 0.5f);

      float span_min = std::max(1e-6f, max_dmin - min_dmin);
      float span_avg = std::max(1e-6f, max_davg - min_davg);

      for (size_t k = 0; k < idx.size(); ++k)
      {
        dist_min[k] += rand_amount * span_min * dist_uni(rng);
        dist_avg[k] += rand_amount * span_avg * dist_uni(rng);
      }
    }

    std::vector<ps::Point<float, 2>> data = ps::merge_by_dimension<float, 2>(
        {dist_min, dist_avg});

    size_t k_clusters = std::min(class_count, actual_count);
    auto [centroids, labels] = ps::kmeans_clustering(data, k_clusters);

    // populate item objects with assigned class tags and default radius
    for (size_t i = 0; i < actual_count; ++i)
    {
      size_t   cluster_id = labels[i];
      uint32_t class_id = static_cast<uint32_t>(cluster_id);

      candidate_items.emplace_back(samples[0][i],
                                   samples[1][i],
                                   0.0f,
                                   class_id,
                                   options.default_radius);
    }
  }

  return ScatterField(std::move(candidate_items));
}

} // namespace hmap
