/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/graph.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/vectors.hpp"

#include <unordered_map>

namespace hmap
{

Graph generate_network_alpha_model(const std::vector<float> &xc,
                                   const std::vector<float> &yc,
                                   const std::vector<float> &size,
                                   glm::vec4                 bbox,
                                   const Array              &z,
                                   std::uint32_t             seed,
                                   float                     alpha,
                                   int                       n_dummy_nodes,
                                   float                     dz_weight,
                                   const Array              *p_weight)
{
  if (!validate_non_empty(z)) return Graph();
  if (xc.size() < 2 || xc.size() != yc.size() || xc.size() != size.size())
    return Graph();
  if (p_weight && !validate_same_shape(z, *p_weight)) return Graph();

  size_t nc = xc.size();

  // --- Tesselation: randomly add "dummy" nodes and use Delaunay triangulation
  // to create a mesh

  Graph graph = Graph();
  {
    const glm::vec2 jitter_amount = {0.5f, 0.5f};
    const glm::vec2 stagger_ratio = {0.f, 0.f};

    Cloud cloud = random_cloud_jittered(n_dummy_nodes,
                                        jitter_amount,
                                        stagger_ratio,
                                        seed,
                                        bbox);

    for (size_t k = 0; k < nc; k++)
    {
      Point p = Point(xc[k], yc[k], size[k]);
      cloud.add_point(p);
    }

    // Delaunay triangulation
    graph = cloud.to_graph_delaunay();
    graph.set_values_from_array(z, bbox);
    graph.update_connectivity();
  }

  // --- Road weights

  std::unordered_map<std::uint64_t, float> is_road;
  is_road.reserve(graph.get_nedges());

  // define number of trips between each cities
  size_t             n_trips = nc * (nc - 1) / 2;
  std::vector<float> ntrips;
  std::vector<int>   trips_istart;
  std::vector<int>   trips_iend;
  ntrips.reserve(n_trips);
  trips_istart.reserve(n_trips);
  trips_iend.reserve(n_trips);

  for (size_t i = 0; i < nc; i++)
    for (size_t j = i + 1; j < nc; j++)
    {
      float dist = (xc[i] - xc[j]) * (xc[i] - xc[j]) +
                   (yc[i] - yc[j]) * (yc[i] - yc[j]);

      ntrips.push_back(size[i] * size[j] / (1.f + dist));
      trips_istart.push_back(static_cast<int>(i));
      trips_iend.push_back(static_cast<int>(j));
    }

  // compute adjacency matrix based on the Euclidian distance
  // between points and add elevation difference
  graph.update_adjacency_matrix();

  std::vector<float> local_weight(graph.size());
  if (p_weight != nullptr)
    local_weight = interpolate_values_from_array(graph, *p_weight, bbox);

  for (size_t i = 0; i < graph.size(); i++)
    for (size_t r = 0; r < graph.connectivity[i].size(); r++)
    {
      int j = graph.connectivity[i][r];
      if (j > static_cast<int>(i))
      {
        float dz = graph.points[i].v - graph.points[j].v;
        graph.adjacency_matrix[{static_cast<int>(i), j}] += std::abs(dz) *
                                                            dz_weight;
        graph.adjacency_matrix[{static_cast<int>(i), j}] += local_weight[i] +
                                                            local_weight[j];
        graph.adjacency_matrix[{j, static_cast<int>(i)}] =
            graph.adjacency_matrix[{static_cast<int>(i), j}];
      }
    }

  // start with the most important connections
  std::vector<size_t> ksort = argsort(ntrips);

  for (size_t k = ntrips.size(); k-- > 0;)
  {
    int i0 = static_cast<int>(graph.size() - nc) + trips_istart[ksort[k]];
    int j0 = static_cast<int>(graph.size() - nc) + trips_iend[ksort[k]];

    // shortest path between the two cities (i0 and j0)
    std::vector<int> path = graph.dijkstra(i0, j0);

    if (path.empty()) continue;

    // update road status and discount edge weight on first activation
    for (size_t i = 0; i < path.size() - 1; i++)
    {
      int           i1 = std::min(path[i], path[i + 1]);
      int           i2 = std::max(path[i], path[i + 1]);
      std::uint64_t key = (static_cast<std::uint64_t>(i1) << 32) |
                          static_cast<std::uint32_t>(i2);
      float &road_count = is_road[key];
      road_count += 1.f;

      if (road_count == 1.f)
      {
        graph.adjacency_matrix[{i1, i2}] *= alpha;
        graph.adjacency_matrix[{i2, i1}] = graph.adjacency_matrix[{i1, i2}];
      }
    }
  }

  // --- Remove orphan edges and rebuild road network graph

  Graph network = Graph(graph.get_x(), graph.get_y());

  for (size_t i = 0; i < graph.size(); i++)
    for (size_t r = 0; r < graph.connectivity[i].size(); r++)
    {
      int j = graph.connectivity[i][r];
      if (j > static_cast<int>(i))
      {
        std::uint64_t key = (static_cast<std::uint64_t>(i) << 32) |
                            static_cast<std::uint32_t>(j);
        auto it = is_road.find(key);
        if (it != is_road.end() && it->second > 0.f)
          network.add_edge({static_cast<int>(i), j}, it->second);
      }
    }

  // store city size in node value (equals to 0 if the node is not a
  // city)
  for (size_t i = 0; i < network.size(); i++)
    if (i < network.size() - nc)
      network.points[i].v = 0.f;
    else
      network.points[i].v = size[i - network.size() + nc];

  // final clean-up
  network = network.remove_orphan_points();
  network.update_adjacency_matrix();

  return network;
}

} // namespace hmap
