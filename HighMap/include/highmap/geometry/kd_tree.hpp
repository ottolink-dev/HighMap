/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file kd_tree.hpp
 * @brief 2D KD-tree for spatial point queries (nearest neighbor and radius
 * search).
 */
#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <glm/glm.hpp>

namespace hmap
{

class Cloud;
struct Point;

/**
 * @brief 2D KD-Tree for spatial point queries.
 *
 * Supports construction from separate coordinate arrays (x, y), point vectors
 * (std::vector<Point>), and Cloud objects.
 *
 * @note Input point data must remain valid for the lifetime of this KDTree
 * object.
 */
class KDTree
{
public:
  // ==========================================================================
  //  Constructors & Destructor
  // ==========================================================================

  /**
   * @brief Construct an empty KD-tree.
   */
  KDTree();

  /**
   * @brief Destructor.
   */
  ~KDTree();

  /**
   * @brief Construct a KD-tree referencing separate coordinate arrays.
   * @param x Vector of x coordinates.
   * @param y Vector of y coordinates.
   */
  KDTree(const std::vector<float> &x, const std::vector<float> &y);

  /**
   * @brief Construct a KD-tree referencing a vector of points.
   * @param points Vector of Point objects.
   */
  explicit KDTree(const std::vector<Point> &points);

  /**
   * @brief Construct a KD-tree referencing a Cloud.
   * @param cloud Cloud object.
   */
  explicit KDTree(const Cloud &cloud);

  /**
   * @brief Move constructor.
   */
  KDTree(KDTree &&other) noexcept;

  /**
   * @brief Move assignment operator.
   */
  KDTree &operator=(KDTree &&other) noexcept;

  KDTree(const KDTree &) = delete;
  KDTree &operator=(const KDTree &) = delete;

  // ==========================================================================
  //  Capacity & Status
  // ==========================================================================

  /**
   * @brief Check whether the KD-tree is empty.
   * @return true if empty, false otherwise.
   */
  bool empty() const noexcept;

  /**
   * @brief Get the number of points indexed by the KD-tree.
   * @return Number of points.
   */
  size_t size() const noexcept;

  // ==========================================================================
  //  Queries
  // ==========================================================================

  /**
   * @brief Estimate the range of neighbor distances across the dataset.
   *
   * Computes the minimum and maximum Euclidean distance to the k-th nearest
   * neighbor across all indexed points.
   *
   * @param  k_neighbors Number of neighbors considered.
   * @return             vec2(min_euclidean_distance, max_euclidean_distance)
   */
  glm::vec2 compute_neighbor_distance_range(size_t k_neighbors) const;

  /**
   * @brief Find the index of the nearest neighbor to a 2D query position.
   *
   * Fast 1-nearest-neighbor query that performs no heap allocations.
   *
   * @param  x_query Query x coordinate.
   * @param  y_query Query y coordinate.
   * @return         Index of the nearest neighbor point (0 if empty).
   */
  size_t nearest(float x_query, float y_query) const;

  /**
   * @brief Find the index of the nearest neighbor to a Point.
   */
  size_t nearest(const Point &p) const;

  /**
   * @brief Find the index of the nearest neighbor to a glm::vec2.
   */
  size_t nearest(const glm::vec2 &xy) const;

  /**
   * @brief Find the index and squared distance of the nearest neighbor.
   *
   * @param  x_query Query x coordinate.
   * @param  y_query Query y coordinate.
   * @return         Pair of (index, squared_distance).
   */
  std::pair<size_t, float> nearest_with_distance_squared(float x_query,
                                                         float y_query) const;

  /**
   * @brief Find the index and squared distance of the nearest neighbor to a
   * Point.
   */
  std::pair<size_t, float> nearest_with_distance_squared(const Point &p) const;

  /**
   * @brief Find the index and squared distance of the nearest neighbor to a
   * glm::vec2.
   */
  std::pair<size_t, float> nearest_with_distance_squared(
      const glm::vec2 &xy) const;

  /**
   * @brief Perform a k-nearest neighbor search.
   *
   * @param      x_query     Query x coordinate.
   * @param      y_query     Query y coordinate.
   * @param      k_neighbors Number of neighbors to retrieve.
   * @param[out] indices     Output indices of neighbors.
   * @param[out] distances   Output squared distances to neighbors.
   */
  void neighbor_search(float                x_query,
                       float                y_query,
                       size_t               k_neighbors,
                       std::vector<size_t> &indices,
                       std::vector<float>  &distances) const;

  /**
   * @brief Perform a k-nearest neighbor search for a Point.
   */
  void neighbor_search(const Point         &p,
                       size_t               k_neighbors,
                       std::vector<size_t> &indices,
                       std::vector<float>  &distances) const;

  /**
   * @brief Perform a k-nearest neighbor search for a glm::vec2.
   */
  void neighbor_search(const glm::vec2     &xy,
                       size_t               k_neighbors,
                       std::vector<size_t> &indices,
                       std::vector<float>  &distances) const;

  /**
   * @brief Perform a radius-based neighbor search.
   *
   * @param  x_query Query x coordinate.
   * @param  y_query Query y coordinate.
   * @param  radius  Search radius (Euclidean distance).
   * @return         Vector of (index, squared_distance) pairs.
   */
  std::vector<std::pair<size_t, float>> radius_search(float x_query,
                                                      float y_query,
                                                      float radius) const;

  /**
   * @brief Perform a radius-based neighbor search for a Point.
   */
  std::vector<std::pair<size_t, float>> radius_search(const Point &p,
                                                      float radius) const;

  /**
   * @brief Perform a radius-based neighbor search for a glm::vec2.
   */
  std::vector<std::pair<size_t, float>> radius_search(const glm::vec2 &xy,
                                                      float radius) const;

  struct Impl;

private:
  std::unique_ptr<Impl> p_impl;
};

} // namespace hmap