/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file forest_seeding.hpp
 * @copyright Copyright (c) 2025 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/flora/forest.hpp"
#include "highmap/flora/species.hpp"

namespace hmap
{

/**
 * @struct ForestSeedingOptions
 * @brief Configuration parameters for forest seeding algorithms.
 */
struct ForestSeedingOptions
{
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}; ///< Bounding box {xmin, xmax,
  // ymin, ymax}.
  uint32_t seed = 0;                   ///< Random number generator seed.
  float    exclusion_threshold = 0.5f; ///< Threshold for exclusion map
  // masking.
  std::vector<Species> species = {}; ///< Optional Species definitions.
};

// ============================================================================
//  Seeding Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Generates a clustered forest distribution using inverse transform
 * sampling for cluster centers and local cluster regions.
 *
 * @param  species_count      Number of distinct tree species.
 * @param  tree_count         Target total number of trees.
 * @param  density            Global density map array.
 * @param  exclusion          Exclusion map array.
 * @param  cluster_spread     Half-extent of the bounding box surrounding each
 *                            cluster center.
 * @param  points_per_cluster Number of child points sampled per parent cluster.
 * @param  options            Seeding options.
 * @return                    Forest Sampled forest container.
 */
Forest seed_forest_clusters(size_t                      species_count,
                            size_t                      tree_count,
                            const Array                &density,
                            const Array                &exclusion,
                            float                       cluster_spread = 0.05f,
                            size_t                      points_per_cluster = 16,
                            const ForestSeedingOptions &options = {});

/**
 * @brief Generates a multi-species forest distribution using 2D inverse
 * transform sampling for positions and PointSampler k-means clustering to tag
 * species.
 *
 * @note Species abundance is governed by spatial cluster partitioning in the
 * compactness feature space; species weights in @p options.species are not
 * taken into account by this function.
 *
 * @param  species_count      Number of distinct tree species.
 * @param  tree_count         Target total number of trees.
 * @param  density            Global density map array.
 * @param  exclusion          Exclusion map array.
 * @param  cluster_randomness Randomness factor in [0, 1] added to cluster
 *                            features.
 * @param  options            Seeding options.
 * @return                    Forest Sampled forest container.
 */
Forest seed_forest_kmeans(size_t                      species_count,
                          size_t                      tree_count,
                          const Array                &density,
                          const Array                &exclusion,
                          float                       cluster_randomness = 0.0f,
                          const ForestSeedingOptions &options = {});

} // namespace hmap
