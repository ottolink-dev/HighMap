/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
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

#include "highmap/flora/forest_seeding.hpp"
#include "highmap/geometry/inverse_sampler_2d.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

// ============================================================================
//  Seeding Functions (Alphabetically Sorted)
// ============================================================================

Forest seed_forest_clusters(size_t                      species_count,
                            size_t                      tree_count,
                            const Array                &density,
                            const Array                &exclusion,
                            float                       cluster_spread,
                            size_t                      points_per_cluster,
                            const ForestSeedingOptions &options)
{
  // validate inputs
  if (species_count == 0 || tree_count == 0 || points_per_cluster == 0 ||
      !validate_non_empty(density))
    return Forest();

  float width = options.bbox.y - options.bbox.x;
  float height = options.bbox.w - options.bbox.z;
  if (width <= 1e-7f || height <= 1e-7f) return Forest();

  // --- Modulate Density by Exclusion Map

  Array effective_density = density;
  bool  has_exclusion = !exclusion.vector.empty() || exclusion.shape.x > 0;

  if (has_exclusion)
  {
    if (!validate_non_empty(exclusion) ||
        !validate_same_shape(density, exclusion))
      return Forest();

    for (int j = 0; j < density.shape.y; ++j)
      for (int i = 0; i < density.shape.x; ++i)
      {
        if (exclusion(i, j) >= options.exclusion_threshold)
          effective_density(i, j) = 0.0f;
      }
  }

  InverseSampler2D sampler(effective_density, options.seed, options.bbox);
  if (sampler.get_total_weight() <= 1e-7f) return Forest();

  // --- Determine Species Quotas

  std::vector<float> effective_weights;
  std::vector<float> effective_radii;

  if (!options.species.empty())
  {
    effective_weights.reserve(options.species.size());
    effective_radii.reserve(options.species.size());
    for (const auto &sp : options.species)
    {
      effective_weights.push_back(sp.weight);
      effective_radii.push_back(sp.radius);
    }
  }

  std::vector<size_t> species_tree_counts(species_count, 0);
  if (effective_weights.size() == species_count)
  {
    float total_w = 0.0f;
    for (float w : effective_weights)
      total_w += std::max(0.0f, w);

    if (total_w > 1e-6f)
    {
      size_t allocated = 0;
      for (size_t s = 0; s < species_count; ++s)
      {
        float norm_w = std::max(0.0f, effective_weights[s]) / total_w;
        species_tree_counts[s] = static_cast<size_t>(
            std::round(static_cast<float>(tree_count) * norm_w));
        allocated += species_tree_counts[s];
      }

      // adjust difference due to rounding
      if (allocated < tree_count)
        species_tree_counts[0] += (tree_count - allocated);
      else if (allocated > tree_count &&
               species_tree_counts[0] >= (allocated - tree_count))
        species_tree_counts[0] -= (allocated - tree_count);
    }
  }
  else
  {
    size_t base_count = tree_count / species_count;
    size_t remainder = tree_count % species_count;
    for (size_t s = 0; s < species_count; ++s)
      species_tree_counts[s] = base_count + (s < remainder ? 1 : 0);
  }

  // --- Generate Clusters per Species

  std::vector<Tree> candidate_trees;
  candidate_trees.reserve(tree_count);

  for (size_t s = 0; s < species_count; ++s)
  {
    if (species_tree_counts[s] == 0) continue;

    size_t target_s = species_tree_counts[s];
    size_t cluster_count = std::max<size_t>(
        1,
        (target_s + points_per_cluster - 1) / points_per_cluster);

    auto parent_samples = sampler.sample(cluster_count);

    float default_radius = (s < effective_radii.size())
                               ? effective_radii[s]
                               : HMAP_DEFAULT_TREE_RADIUS;

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
        candidate_trees.emplace_back(pt.x,
                                     pt.y,
                                     0.0f,
                                     static_cast<uint32_t>(s),
                                     default_radius);
        count_for_s++;
      }
    }
  }

  Forest forest(std::move(candidate_trees));
  return forest;
}

Forest seed_forest_kmeans(size_t                      species_count,
                          size_t                      tree_count,
                          const Array                &density,
                          const Array                &exclusion,
                          float                       cluster_randomness,
                          const ForestSeedingOptions &options)
{
  // validate inputs
  if (species_count == 0 || tree_count == 0 || !validate_non_empty(density))
    return Forest();

  float width = options.bbox.y - options.bbox.x;
  float height = options.bbox.w - options.bbox.z;
  if (width <= 1e-7f || height <= 1e-7f) return Forest();

  // --- Modulate Density by Exclusion Map

  Array effective_density = density;
  bool  has_exclusion = !exclusion.vector.empty() || exclusion.shape.x > 0;

  if (has_exclusion)
  {
    if (!validate_non_empty(exclusion) ||
        !validate_same_shape(density, exclusion))
      return Forest();

    for (int j = 0; j < density.shape.y; ++j)
      for (int i = 0; i < density.shape.x; ++i)
      {
        if (exclusion(i, j) >= options.exclusion_threshold)
          effective_density(i, j) = 0.0f;
      }
  }

  InverseSampler2D sampler(effective_density, options.seed, options.bbox);
  if (sampler.get_total_weight() <= 1e-7f) return Forest();

  // --- Sample Coordinates Using 2D Inverse Sampling

  std::vector<float> effective_radii;
  if (!options.species.empty())
  {
    effective_radii.reserve(options.species.size());
    for (const auto &sp : options.species)
      effective_radii.push_back(sp.radius);
  }

  // generate candidate point positions based on the 2d probability density
  auto   samples = sampler.sample(tree_count);
  size_t actual_count = samples[0].size();
  if (actual_count == 0) return Forest();

  // --- Cluster and Tag Species

  std::vector<Tree> candidate_trees;
  candidate_trees.reserve(actual_count);

  if (species_count == 1 || actual_count == 1)
  {
    // single species: assign all sampled points to species 0
    float default_radius = (!effective_radii.empty())
                               ? effective_radii[0]
                               : HMAP_DEFAULT_TREE_RADIUS;

    for (size_t i = 0; i < actual_count; ++i)
    {
      candidate_trees.emplace_back(samples[0][i],
                                   samples[1][i],
                                   0.0f,
                                   0u,
                                   default_radius);
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

    // extract local compactness features using k-nearest neighbors:
    // rather than clustering on raw (x, y) coordinates (which would partition
    // space into geometric voronoi cells), we compute the minimum distance
    // (dmin) and average distance (davg) to nearest neighbors. this
    // characterizes whether a point is in a dense cluster core, transition
    // zone, or isolated. note that species abundance emerges from cluster
    // geometry and options.species_weights is not taken into account.
    size_t k_neighbors = 4;
    size_t safe_k_neighbors = std::min(k_neighbors, actual_count - 1);

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

    // merge [dmin, davg] into 2d feature points for k-means partitioning
    std::vector<ps::Point<float, 2>> data = ps::merge_by_dimension<float, 2>(
        {dist_min, dist_avg});

    // cluster points in the compactness feature space (data is normalized
    // internally)
    size_t k_clusters = std::min(species_count, actual_count);
    auto [centroids, labels] = ps::kmeans_clustering(data, k_clusters);

    // establish correspondence between cluster compactness and species radii:
    // clusters with smaller neighbor distances (dmin/davg, dense areas)
    // correspond to species with smaller radii, while clusters with larger
    // neighbor distances (sparse areas) correspond to species with larger
    // radii.
    std::vector<uint32_t> cluster_to_species(k_clusters);
    for (size_t c = 0; c < k_clusters; ++c)
    {
      cluster_to_species[c] = static_cast<uint32_t>(c);
    }

    if (!effective_radii.empty() && effective_radii.size() == k_clusters)
    {
      // sort cluster indices by centroid distance in compactness feature space
      std::vector<size_t> cluster_order(k_clusters);
      std::iota(cluster_order.begin(), cluster_order.end(), 0);
      std::sort(cluster_order.begin(),
                cluster_order.end(),
                [&centroids](size_t a, size_t b)
                {
                  float norm_a = centroids[a][0] + centroids[a][1];
                  float norm_b = centroids[b][0] + centroids[b][1];
                  return norm_a < norm_b;
                });

      // sort available species indices by their radii ascending
      size_t num_species = std::min(species_count, effective_radii.size());
      std::vector<size_t> species_order(num_species);
      std::iota(species_order.begin(), species_order.end(), 0);
      std::sort(species_order.begin(),
                species_order.end(),
                [&effective_radii](size_t a, size_t b)
                { return effective_radii[a] < effective_radii[b]; });

      for (size_t rank = 0; rank < k_clusters; ++rank)
      {
        size_t c_idx = cluster_order[rank];
        size_t sp_idx = (rank < num_species) ? species_order[rank] : rank;
        cluster_to_species[c_idx] = static_cast<uint32_t>(sp_idx);
      }
    }

    // populate tree objects with assigned species tags and radii
    for (size_t i = 0; i < actual_count; ++i)
    {
      size_t   cluster_id = labels[i];
      uint32_t species_id = (cluster_id < cluster_to_species.size())
                                ? cluster_to_species[cluster_id]
                                : static_cast<uint32_t>(cluster_id);
      float    radius = (species_id < effective_radii.size())
                            ? effective_radii[species_id]
                            : HMAP_DEFAULT_TREE_RADIUS;

      candidate_trees.emplace_back(samples[0][i],
                                   samples[1][i],
                                   0.0f,
                                   species_id,
                                   radius);
    }
  }

  Forest forest(std::move(candidate_trees));
  return forest;
}

} // namespace hmap
