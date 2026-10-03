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
#include <optional>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/flora/tree.hpp"
#include "highmap/geometry/cloud.hpp"

namespace hmap
{

struct Species;

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
         float        default_radius = HMAP_DEFAULT_TREE_RADIUS);

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
   * @brief Densifies the forest by adding Voronoi vertices (circumcenters of
   * Delaunay triangles).
   *
   * Computes the Delaunay triangulation of the existing trees. For each
   * triangle, computes its circumcenter (a Voronoi vertex) and inserts a new
   * tree with the majority species identifier among the triangle's 3 vertices
   * and the specified default radius.
   *
   * @param default_radius Default radius assigned to newly added trees.
   */
  void densify(float default_radius = HMAP_DEFAULT_TREE_RADIUS);

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
   * @brief Prunes trees in the forest using density-based acceptance sampling
   * modulated to reach an approximate target retention ratio.
   *
   * @note The retention ratio is achieved in expectation (approximate) via
   * probabilistic Bernoulli sampling biased by local density values.
   *
   * @param density_mask 2D array defining the spatial density field.
   * @param target_ratio Approximate target fraction of trees to retain in [0,
   * 1].
   * @param seed         Random seed for reproducible pruning.
   * @param bbox         Bounding box defining the domain of the density mask.
   */
  void prune_density(const Array     &density_mask,
                     float            target_ratio = 0.8f,
                     uint32_t         seed = 0,
                     const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Prunes unviable trees whose distance to their nearest neighbor is
   * insufficient to meet their species' minimum radius requirement or pairwise
   * collision distance.
   *
   * Uses a KD-Tree / nearest-neighbor search to find the distance $d_{nn}$ to
   * the closest neighbor. If estimated radius $\alpha \cdot d_{nn} < r_{\min}$
   * (or if overlapping with a larger tree), the tree is pruned.
   *
   * @param species          Vector of species definitions (with radius_min,
   *                         radius_max, competition_factor).
   * @param prune_collisions If true, also resolves overlapping tree crowns by
   *                         discarding the smaller tree.
   */
  void prune_unviable(const std::vector<Species> &species = {},
                      bool                        prune_collisions = true);

  /**
   * @brief Reinforces spatial clustering of species by iteratively assigning
   * each tree the dominant species among its nearest neighbors.
   *
   * @param iterations   Number of smoothing/reinforcement iterations.
   * @param k_neighbors  Number of spatial nearest neighbors to query.
   * @param include_self If true, considers the tree's own current species in
   * the majority vote.
   */
  void reinforce_species_clusters(size_t iterations = 2,
                                  size_t k_neighbors = 4,
                                  bool   include_self = true);

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
   * @brief Randomly shuffles species identifiers between neighboring trees.
   *
   * For a randomly selected ratio of trees, a tree swaps its species with one
   * of its @p k_neighbors. If all neighbors have the same species as the
   * candidate node, no swap is made.
   *
   * @param ratio       Fraction of trees to attempt species shuffling on in [0,
   *                    1].
   * @param k_neighbors Number of nearest spatial neighbors to consider.
   * @param seed        Random seed for reproducibility.
   */
  void shuffle_species(float    ratio = 0.1f,
                       size_t   k_neighbors = 4,
                       uint32_t seed = 0);

  /**
   * @brief Converts the forest into a Cloud of 2D points (x, y) with value set
   * to radius.
   * @return Cloud Point cloud representation.
   */
  Cloud to_cloud() const;

  /**
   * @brief Exports the local density map of the forest as a 2D Array using 2D
   * Kernel Density Estimation (Gaussian splatting / KDE).
   *
   * Splats each tree into the regular grid and accumulates smooth 2D Gaussian
   * kernels with bandwidth @p sigma. If @p weighted_by_crown is true, each
   * contribution is weighted by the tree crown surface ($\pi r^2$) instead of
   * unit count.
   *
   * @param  shape             Dimensions of the output density map {nx, ny}.
   * @param  sigma             Gaussian kernel standard deviation (bandwidth) in
   *                           world/domain coordinates. If <= 0, automatically
   *                           derived from the bounding box and grid shape.
   * @param  species_id        Optional species identifier to compute density
   * for a specific species only.
   * @param  weighted_by_crown If true, weights tree density by crown surface
   *                           ($\pi r^2$).
   * @param  bbox              Bounding box defining the domain {xmin, xmax,
   * ymin, ymax}.
   * @return                   Array 2D continuous density map array.
   */
  Array to_density_map(glm::ivec2              shape,
                       float                   sigma = 0.0f,
                       std::optional<uint32_t> species_id = std::nullopt,
                       bool                    weighted_by_crown = false,
                       glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}) const;

  /**
   * @brief Exports the forest trees to a CSV file (x, y, z, species_id,
   * radius).
   * @param fname Output file path.
   */
  void to_csv(const std::string &fname) const;

  /**
   * @brief Exports a visual representation of the forest as a PNG image.
   *
   * @param fname      Output PNG file path.
   * @param shape      Image dimensions {width, height}.
   * @param background Optional background terrain array (rendered as grayscale,
   *                   resampled if needed).
   * @param bbox       Bounding box {xmin, xmax, ymin, ymax} of the domain.
   */
  void to_png(const std::string &fname,
              glm::ivec2         shape,
              const Array       &background = {},
              glm::vec4          bbox = {0.f, 1.f, 0.f, 1.f}) const;

  /**
   * @brief Returns a multi-line formatted summary string of the forest.
   * @return std::string Pretty-printed summary.
   */
  std::string to_string() const;

protected:
  std::vector<Tree> trees = {}; ///< List of tree instances.
};

} // namespace hmap
