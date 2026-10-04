/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file scatter_seeding.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/scatter/scatter_field.hpp"

namespace hmap
{

/**
 * @struct ScatterSeedingOptions
 * @brief Configuration parameters for scatter seeding algorithms.
 */
struct ScatterSeedingOptions
{
  glm::vec4 bbox = {0.f,
                    1.f,
                    0.f,
                    1.f}; ///< Bounding box {xmin, xmax, ymin, ymax}.
  uint32_t  seed = 0;     ///< Random number generator seed.
  float     exclusion_threshold = 0.5f; ///< Threshold for exclusion map masking.
  std::vector<float> class_weights =
      {}; ///< Relative abundance weights per class_id.
  float default_radius = HMAP_DEFAULT_SCATTER_RADIUS; ///< Default radius.
};

// ============================================================================
//  Seeding Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Generates a clustered scatter distribution using inverse transform
 * sampling for cluster centers and local cluster regions.
 *
 * @param  class_count        Number of distinct classes/categories.
 * @param  item_count         Target total number of items.
 * @param  density            Global density map array.
 * @param  exclusion          Exclusion map array.
 * @param  cluster_spread     Half-extent of the bounding box surrounding each
 *                            cluster center.
 * @param  points_per_cluster Number of child points sampled per parent cluster.
 * @param  options            Seeding options.
 * @return                    ScatterField Sampled scatter field container.
 */
ScatterField seed_scatter_clusters(size_t       class_count,
                                   size_t       item_count,
                                   const Array &density,
                                   const Array &exclusion,
                                   float        cluster_spread = 0.05f,
                                   size_t       points_per_cluster = 8,
                                   const ScatterSeedingOptions &options = {});

/**
 * @brief Generates a multi-class scatter distribution using 2D inverse
 * transform sampling for positions and PointSampler k-means clustering to tag
 * classes.
 *
 * @param  class_count        Number of distinct item classes.
 * @param  item_count         Target total number of items.
 * @param  density            Global density map array.
 * @param  exclusion          Exclusion map array.
 * @param  cluster_randomness Randomness factor in [0, 1] added to cluster
 *                            features.
 * @param  k_neighbors        Number of nearest spatial neighbors used to
 *                            extract compactness features.
 * @param  options            Seeding options.
 * @return                    ScatterField Sampled scatter field container.
 */
ScatterField seed_scatter_kmeans(size_t       class_count,
                                 size_t       item_count,
                                 const Array &density,
                                 const Array &exclusion,
                                 float        cluster_randomness = 0.0f,
                                 size_t       k_neighbors = 4,
                                 const ScatterSeedingOptions &options = {});

} // namespace hmap
