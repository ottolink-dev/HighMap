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
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

/**
 * @struct ForestSeedingOptions
 * @brief Configuration parameters for forest seeding algorithms.
 */
struct ForestSeedingOptions
{
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}; ///< Bounding box {xmin, xmax, ymin, ymax}.
  uint32_t seed = 0;                     ///< Random number generator seed.
  float    exclusion_threshold = 0.5f;   ///< Threshold for exclusion map masking.
  std::vector<Species> species = {};     ///< Optional Species definitions.
};

// ============================================================================
//  Seeding Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Computes an environmental suitability criteria map in [0, 1]
 * representing tree probability density based on elevation, slope talus, aspect
 * angle, topographic wetness index (TWI), an optional secondary density map,
 * and an optional exclusion mask.
 *
 * @param  heightmap             Input terrain heightmap.
 * @param  min_elev              Minimum elevation threshold.
 * @param  max_elev              Maximum elevation threshold.
 * @param  elev_transition_width Elevation transition smoothing width.
 * @param  min_talus             Slope/talus threshold for full suitability
 *                               (1.0).
 * @param  max_talus             Slope/talus threshold for zero suitability
 *                               (0.0).
 * @param  angle                 Preferred slope aspect angle in degrees.
 * @param  angle_width           Slope aspect angle tolerance half-width in
 *                               degrees.
 * @param  weight_elev           Weight for elevation criterion in [0, 1].
 * @param  weight_talus          Weight for talus/slope criterion in [0, 1].
 * @param  weight_angle          Weight for aspect angle criterion in [0, 1].
 * @param  weight_twi            Weight for topographic wetness index criterion
 *                               in [0, 1].
 * @param  secondary_density     Optional secondary density or noise modulation
 *                               map.
 * @param  weight_secondary      Weight for optional secondary density map in
 *                               [0, 1].
 * @param  exclusion_mask        Optional exclusion mask (where <= 0, density is
 *                               set to 0).
 * @return                       Array Computed density array in [0, 1].
 */
Array build_tree_density(const Array &heightmap,
                         float        min_elev = 0.0f,
                         float        max_elev = 0.9f,
                         float        elev_transition_width = 0.1f,
                         float        min_talus = 0.003f,
                         float        max_talus = 0.005f,
                         float        angle = 30.0f,
                         float        angle_width = 90.0f,
                         float        weight_elev = 1.0f,
                         float        weight_talus = 1.0f,
                         float        weight_angle = 0.8f,
                         float        weight_twi = 0.5f,
                         const Array &secondary_density = {},
                         float        weight_secondary = 1.0f,
                         const Array &exclusion_mask = {});

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
                            size_t                      points_per_cluster = 8,
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
 * @param  k_neighbors        Number of nearest spatial neighbors used to
 *                            extract compactness features.
 * @param  options            Seeding options.
 * @return                    Forest Sampled forest container.
 */
Forest seed_forest_kmeans(size_t                      species_count,
                          size_t                      tree_count,
                          const Array                &density,
                          const Array                &exclusion,
                          float                       cluster_randomness = 0.0f,
                          size_t                      k_neighbors = 4,
                          const ForestSeedingOptions &options = {});

} // namespace hmap

namespace hmap::va
{

/**
 * @brief Generates a clustered forest distribution across a tiled VirtualArray
 * domain using inverse transform sampling and bounding-box-aware tile merging.
 *
 * @param  species_count      Number of distinct tree species.
 * @param  tree_count         Target total number of trees.
 * @param  density            Global density VirtualArray.
 * @param  exclusion          Optional exclusion map VirtualArray.
 * @param  cluster_spread     Half-extent of the bounding box surrounding each
 *                            cluster center.
 * @param  points_per_cluster Number of child points sampled per parent cluster.
 * @param  options            Seeding options.
 * @param  cm                 Virtual array compute mode.
 * @return                    Forest Sampled forest container.
 */
Forest seed_forest_clusters(size_t                      species_count,
                            size_t                      tree_count,
                            const VirtualArray         &density,
                            const VirtualArray         &exclusion = {},
                            float                       cluster_spread = 0.05f,
                            size_t                      points_per_cluster = 8,
                            const ForestSeedingOptions &options = {},
                            const ComputeMode          &cm = {});

} // namespace hmap::va
