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

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/flora/tree.hpp"
#include "highmap/geometry/cloud.hpp"

namespace hmap
{

/**
 * @class Forest
 * @brief Container class representing a collection of Tree instances.
 */
class Forest
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
   * @brief Constructs a Forest from a Cloud of 2D points.
   *
   * If a point's value `v` is non-zero, it is used as the tree's radius;
   * otherwise @p default_radius is assigned.
   *
   * @param cloud          Input point cloud.
   * @param species_id     Species identifier assigned to all imported trees.
   * @param default_radius Default radius assigned if point value is 0.
   */
  Forest(const Cloud &cloud,
         uint32_t     species_id = 0,
         float        default_radius = 1.0f);

  // ==========================================================================
  //  Container Interface
  // ==========================================================================

  Tree       &at(size_t index);
  const Tree &at(size_t index) const;

  Tree       &back();
  const Tree &back() const;

  auto begin() noexcept
  {
    return trees.begin();
  }
  auto begin() const noexcept
  {
    return trees.begin();
  }
  auto cbegin() const noexcept
  {
    return trees.cbegin();
  }

  size_t capacity() const noexcept
  {
    return trees.capacity();
  }
  void clear() noexcept
  {
    trees.clear();
  }

  Tree *data() noexcept
  {
    return trees.data();
  }
  const Tree *data() const noexcept
  {
    return trees.data();
  }

  template <typename... Args> Tree &emplace_back(Args &&...args)
  {
    return trees.emplace_back(std::forward<Args>(args)...);
  }

  bool empty() const noexcept
  {
    return trees.empty();
  }

  auto end() noexcept
  {
    return trees.end();
  }
  auto end() const noexcept
  {
    return trees.end();
  }
  auto cend() const noexcept
  {
    return trees.cend();
  }

  Tree       &front();
  const Tree &front() const;

  Tree &operator[](size_t index)
  {
    return trees[index];
  }
  const Tree &operator[](size_t index) const
  {
    return trees[index];
  }

  void push_back(const Tree &tree)
  {
    trees.push_back(tree);
  }
  void push_back(Tree &&tree)
  {
    trees.push_back(std::move(tree));
  }

  void reserve(size_t new_cap)
  {
    trees.reserve(new_cap);
  }
  size_t size() const noexcept
  {
    return trees.size();
  }

  // ==========================================================================
  //  Operations
  // ==========================================================================

  /**
   * @brief Returns a new Forest containing only trees matching the given
   * species.
   * @param  species_id Species identifier to filter by.
   * @return            Forest Filtered forest.
   */
  Forest filter_by_species(uint32_t species_id) const;

  /**
   * @brief Computes the 2D bounding box enclosing all tree positions.
   * @return glm::vec4 Bounding box as {xmin, xmax, ymin, ymax}.
   */
  glm::vec4 get_bbox() const;

  /**
   * @brief Returns a sorted list of unique species identifiers present in the
   * forest.
   * @return std::vector<uint32_t> Unique species IDs.
   */
  std::vector<uint32_t> get_species_ids() const;

  /**
   * @brief Sets the elevation (z-coordinate) of all trees by sampling a terrain
   * heightmap.
   *
   * @param elevation Terrain heightmap array.
   * @param bbox      Bounding box defining the domain of the elevation array.
   */
  void set_elevation_from_terrain(const Array     &elevation,
                                  const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Converts the forest into a Cloud of 2D points (x, y) with value set
   * to radius.
   * @return Cloud Point cloud representation.
   */
  Cloud to_cloud() const;

protected:
  std::vector<Tree> trees = {}; ///< List of tree instances.
};

} // namespace hmap
