/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "highmap/terrain_tri_mesh.hpp"

namespace hmap
{

TerrainTriMesh::ShortestPathResult TerrainTriMesh::compute_shortest_paths(
    size_t start,
    bool   use_delta_z,
    float  elevation_weight) const
{
  constexpr size_t INVALID = std::numeric_limits<size_t>::max();

  ShortestPathResult result;

  const size_t n = points.size();

  result.distance.assign(n, std::numeric_limits<float>::infinity());
  result.parent.assign(n, INVALID);

  using QueueNode = std::pair<float, size_t>;

  std::
      priority_queue<QueueNode, std::vector<QueueNode>, std::greater<QueueNode>>
          queue;

  result.distance[start] = 0.f;
  result.parent[start] = start;

  queue.emplace(0.f, start);

  while (!queue.empty())
  {
    auto [dist, v] = queue.top();
    queue.pop();

    if (dist > result.distance[v]) continue;

    for (const Neighbor &n : neighbors.adjacency[v])
    {
      float edge_cost;

      if (use_delta_z)
      {
        float dz = points[v].z - points[n.index].z;
        edge_cost = std::abs(dz) + elevation_weight * points[v].z;
      }
      else
      {
        edge_cost = n.distance2d;
      }

      float new_dist = dist + edge_cost;

      if (new_dist < result.distance[n.index])
      {
        result.distance[n.index] = new_dist;
        result.parent[n.index] = v;

        queue.emplace(new_dist, n.index);
      }
    }
  }

  return result;
}

TerrainTriMesh::ShortestPathResult TerrainTriMesh::
    compute_shortest_paths_to_hull(bool  use_delta_z,
                                   float elevation_weight) const
{
  constexpr size_t INVALID = std::numeric_limits<size_t>::max();

  ShortestPathResult result;

  const size_t n = points.size();

  result.distance.assign(n, std::numeric_limits<float>::infinity());
  result.parent.assign(n, INVALID);

  using QueueNode = std::pair<float, size_t>;

  std::
      priority_queue<QueueNode, std::vector<QueueNode>, std::greater<QueueNode>>
          queue;

  for (size_t v : convex_hull)
  {
    result.distance[v] = 0.f;
    result.parent[v] = v;
    queue.emplace(0.f, v);
  }

  while (!queue.empty())
  {
    auto [dist, v] = queue.top();
    queue.pop();

    if (dist > result.distance[v]) continue;

    for (const Neighbor &n : neighbors.adjacency[v])
    {
      float edge_cost;

      if (use_delta_z)
      {
        float dz = points[v].z - points[n.index].z;
        edge_cost = std::abs(dz) + elevation_weight * points[v].z;
      }
      else
      {
        edge_cost = n.distance2d;
      }

      float new_dist = dist + edge_cost;

      if (new_dist < result.distance[n.index])
      {
        result.distance[n.index] = new_dist;
        result.parent[n.index] = v;

        queue.emplace(new_dist, n.index);
      }
    }
  }

  return result;
}

std::vector<size_t> TerrainTriMesh::shortest_path(size_t start,
                                                  size_t end,
                                                  bool   use_delta_z,
                                                  float  elevation_weight) const
{
  constexpr size_t INVALID = std::numeric_limits<size_t>::max();

  auto result = compute_shortest_paths(start, use_delta_z, elevation_weight);

  std::vector<size_t> path;

  if (result.parent[end] == INVALID) return path;

  for (size_t v = end;; v = result.parent[v])
  {
    path.push_back(v);

    if (v == start) break;
  }

  std::reverse(path.begin(), path.end());

  return path;
}

std::vector<size_t> TerrainTriMesh::path_to_hull(
    size_t                    start,
    const ShortestPathResult &result) const
{
  constexpr size_t INVALID = std::numeric_limits<size_t>::max();

  std::vector<size_t> path;

  if (result.parent[start] == INVALID) return path;

  for (size_t v = start;; v = result.parent[v])
  {
    path.push_back(v);

    if (result.parent[v] == v) break;
  }

  return path;
}

} // namespace hmap
