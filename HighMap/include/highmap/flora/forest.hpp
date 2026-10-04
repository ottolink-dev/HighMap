/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file forest.hpp
 * @copyright Copyright (c) 2025 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "highmap/flora/species.hpp"
#include "highmap/flora/tree.hpp"
#include "highmap/scatter/scatter_field.hpp"

namespace hmap
{

/**
 * @class Forest
 * @brief Container class representing a collection of Tree instances.
 */
class Forest : public ScatterField
{
public:
  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Default constructor initializing an empty forest.
   */
  Forest() = default;

  /**
   * @brief Constructs a Forest from a vector of trees.
   * @param trees Vector of Tree instances.
   */
  Forest(const std::vector<Tree> &trees);

  /**
   * @brief Move-constructs a Forest from a vector of trees.
   * @param trees Rvalue vector of Tree instances.
   */
  Forest(std::vector<Tree> &&trees) noexcept;

  /**
   * @brief Constructs a Forest from a base ScatterField.
   * @param field Base scatter field.
   */
  Forest(const ScatterField &field);

  /**
   * @brief Move-constructs a Forest from a base ScatterField.
   * @param field Rvalue base scatter field.
   */
  Forest(ScatterField &&field) noexcept;

  /**
   * @brief Constructs a Forest from a Cloud of 2D points.
   *
   * @param cloud          Input point cloud.
   * @param species_id     Species identifier assigned to all imported trees.
   * @param default_radius Default radius assigned if point value is 0.
   */
  Forest(const Cloud &cloud,
         uint32_t     species_id = 0,
         float        default_radius = HMAP_DEFAULT_TREE_RADIUS);

  /**
   * @brief Prunes unviable trees whose distance to their nearest neighbor is
   * insufficient to meet their species' minimum radius requirement or pairwise
   * collision distance.
   *
   * @param species          Vector of species definitions (with radius_min,
   *                         radius_max, competition_factor).
   * @param prune_collisions If true, also resolves overlapping tree crowns by
   *                         discarding the smaller tree.
   */
  void prune_unviable(const std::vector<Species> &species = {},
                      bool                        prune_collisions = true);

  /**
   * @brief Returns a multi-line formatted summary string of the forest.
   * @return std::string Pretty-printed summary.
   */
  std::string to_string() const override;
};

} // namespace hmap
