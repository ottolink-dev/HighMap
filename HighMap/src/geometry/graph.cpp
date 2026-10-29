/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/graph.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"
#include "highmap/sdf.hpp"

#include <format>

namespace hmap
{

const std::vector<Neighbor> Graph::empty_neighbors = {};

// --- Constructors

Graph::Graph(Cloud cloud) : Cloud(std::move(cloud))
{
  this->ensure_adj_capacity(this->size());
}

Graph::Graph(std::vector<Point> points) : Cloud(std::move(points))
{
  this->ensure_adj_capacity(this->size());
}

Graph::Graph(std::vector<float> x, std::vector<float> y)
    : Cloud(std::move(x), std::move(y))
{
  this->ensure_adj_capacity(this->size());
}

// --- Private helpers

void Graph::ensure_adj_capacity(size_t node_count)
{
  if (this->adj_list.size() < node_count) this->adj_list.resize(node_count);
}

// --- Public Graph methods (alphabetical)

void Graph::add_edge(int u, int v, float weight)
{
  size_t max_node = static_cast<size_t>(std::max(u, v) + 1);
  this->ensure_adj_capacity(std::max(this->size(), max_node));

  int edge_idx = static_cast<int>(this->edge_list.size());
  this->edge_list.push_back({u, v, weight});

  if (u >= 0 && static_cast<size_t>(u) < this->adj_list.size())
    this->adj_list[u].push_back({v, weight, edge_idx});
  if (v >= 0 && static_cast<size_t>(v) < this->adj_list.size())
    this->adj_list[v].push_back({u, weight, edge_idx});
}

void Graph::add_edge(int u, int v)
{
  float w = 0.f;
  if (u >= 0 && u < static_cast<int>(this->size()) && v >= 0 &&
      v < static_cast<int>(this->size()))
  {
    w = distance(this->points[u], this->points[v]);
  }
  this->add_edge(u, v, w);
}

void Graph::add_edge(glm::ivec2 edge, float weight)
{
  this->add_edge(edge.x, edge.y, weight);
}

void Graph::add_edge(glm::ivec2 edge)
{
  this->add_edge(edge.x, edge.y);
}

void Graph::clear_edges() noexcept
{
  this->edge_list.clear();
  this->adj_list.clear();
}

std::vector<float> Graph::compute_distance_to_outlet() const
{
  const size_t       n = this->size();
  std::vector<float> dist(n, 0.f);

  if (n == 0) return dist;

  // build directed downstream graph: edge u -> v if elevation(u) > elevation(v)
  // or tie-broken by node index
  std::vector<std::vector<int>> downstream_neighbors(n);
  std::vector<int>              in_degree(n, 0);

  for (const auto &edge : this->edge_list)
  {
    int u = edge.u;
    int v = edge.v;

    if (u < 0 || static_cast<size_t>(u) >= n || v < 0 ||
        static_cast<size_t>(v) >= n || u == v)
      continue;

    float vu = this->points[u].v;
    float vv = this->points[v].v;

    int from = u;
    int to = v;

    if (vu < vv || (vu == vv && u > v))
    {
      from = v;
      to = u;
    }

    downstream_neighbors[from].push_back(to);
    in_degree[to]++;
  }

  // Kahn's algorithm: topological sort from upstream headwaters to downstream
  // outlets
  std::queue<int>  q;
  std::vector<int> in_deg_remaining = in_degree;

  for (size_t i = 0; i < n; ++i)
  {
    if (in_deg_remaining[i] == 0) q.push(static_cast<int>(i));
  }

  std::vector<int> topo_order;
  topo_order.reserve(n);

  while (!q.empty())
  {
    int u = q.front();
    q.pop();
    topo_order.push_back(u);

    for (int v : downstream_neighbors[u])
    {
      in_deg_remaining[v]--;
      if (in_deg_remaining[v] == 0) q.push(v);
    }
  }

  // if any cycle existed (e.g. flat elevations), append remaining nodes
  if (topo_order.size() < n)
  {
    std::vector<bool> visited(n, false);
    for (int u : topo_order)
      visited[u] = true;

    for (size_t i = 0; i < n; ++i)
      if (!visited[i]) topo_order.push_back(static_cast<int>(i));
  }

  // propagate distances to outlet in reverse topological order (from outlets to
  // headwaters)
  for (auto it = topo_order.rbegin(); it != topo_order.rend(); ++it)
  {
    int u = *it;

    if (downstream_neighbors[u].empty())
    {
      dist[u] = 0.f;
      continue;
    }

    float min_d = std::numeric_limits<float>::max();

    for (int v : downstream_neighbors[u])
    {
      float edge_len = glm::distance(
          glm::vec2(this->points[u].x, this->points[u].y),
          glm::vec2(this->points[v].x, this->points[v].y));
      float candidate_d = dist[v] + edge_len;

      if (candidate_d < min_d)
      {
        min_d = candidate_d;
      }
    }

    dist[u] = (min_d < std::numeric_limits<float>::max()) ? min_d : 0.f;
  }

  return dist;
}

std::vector<float> Graph::compute_flow_accumulation() const
{
  const size_t       n = this->size();
  std::vector<float> accum(n, 0.f);

  if (n == 0 || this->edge_list.empty()) return accum;

  // build directed downstream edges: from higher elevation to lower elevation
  std::vector<std::vector<int>> downstream_neighbors(n);
  std::vector<int>              in_degree(n, 0);

  for (const auto &edge : this->edge_list)
  {
    int u = edge.u;
    int v = edge.v;

    if (u < 0 || static_cast<size_t>(u) >= n || v < 0 ||
        static_cast<size_t>(v) >= n || u == v)
      continue;

    float vu = this->points[u].v;
    float vv = this->points[v].v;

    int from = u;
    int to = v;

    if (vu < vv || (vu == vv && u > v))
    {
      from = v;
      to = u;
    }

    downstream_neighbors[from].push_back(to);
    in_degree[to]++;
  }

  // Kahn's algorithm: topological sort from upstream headwaters to downstream
  // outlets
  std::queue<int>  q;
  std::vector<int> in_deg_remaining = in_degree;

  for (size_t i = 0; i < n; ++i)
  {
    if (in_deg_remaining[i] == 0) q.push(static_cast<int>(i));
  }

  std::vector<int> topo_order;
  topo_order.reserve(n);

  while (!q.empty())
  {
    int u = q.front();
    q.pop();
    topo_order.push_back(u);

    for (int v : downstream_neighbors[u])
    {
      in_deg_remaining[v]--;
      if (in_deg_remaining[v] == 0) q.push(v);
    }
  }

  // if any cycle existed, append remaining nodes
  if (topo_order.size() < n)
  {
    std::vector<bool> visited(n, false);
    for (int u : topo_order)
      visited[u] = true;

    for (size_t i = 0; i < n; ++i)
      if (!visited[i]) topo_order.push_back(static_cast<int>(i));
  }

  // accumulate flow in topological order from upstream to downstream
  for (int u : topo_order)
  {
    for (int v : downstream_neighbors[u])
    {
      float edge_len = glm::distance(
          glm::vec2(this->points[u].x, this->points[u].y),
          glm::vec2(this->points[v].x, this->points[v].y));

      accum[v] += (accum[u] + edge_len);
    }
  }

  return accum;
}

std::vector<float> Graph::compute_relative_widths(float exponent) const
{
  const size_t       n = this->size();
  std::vector<float> rel_widths(n, 0.f);

  if (n == 0 || this->edge_list.empty()) return rel_widths;

  // compute downstream distance to outlet for all nodes
  std::vector<float> dist_to_outlet = this->compute_distance_to_outlet();

  // extract directed downstream edges and identify root nodes (in_degree == 0)
  std::vector<std::vector<int>> downstream_adj(n);
  std::vector<int>              in_degree(n, 0);

  for (const auto &edge : this->edge_list)
  {
    int u = edge.u;
    int v = edge.v;

    if (u < 0 || static_cast<size_t>(u) >= n || v < 0 ||
        static_cast<size_t>(v) >= n || u == v)
      continue;

    float vu = this->points[u].v;
    float vv = this->points[v].v;

    int from = u;
    int to = v;

    if (vu < vv || (vu == vv && u > v))
    {
      from = v;
      to = u;
    }

    downstream_adj[from].push_back(to);
    in_degree[to]++;
  }

  // identify the outlet for each node and maximum distance to outlet per
  // drainage tree
  std::vector<int>   outlet_of(n, -1);
  std::vector<float> max_tree_dist(n, 0.f);

  for (size_t i = 0; i < n; ++i)
  {
    size_t curr = i;
    while (!downstream_adj[curr].empty())
    {
      float  min_d = std::numeric_limits<float>::max();
      size_t best_next = curr;

      for (int nxt : downstream_adj[curr])
      {
        if (dist_to_outlet[nxt] < min_d)
        {
          min_d = dist_to_outlet[nxt];
          best_next = static_cast<size_t>(nxt);
        }
      }

      if (best_next == curr) break;
      curr = best_next;
    }

    outlet_of[i] = static_cast<int>(curr);
    max_tree_dist[curr] = std::max(max_tree_dist[curr], dist_to_outlet[i]);
  }

  // compute relative width ratio
  for (size_t i = 0; i < n; ++i)
  {
    // headwater root nodes with no incoming flow start with zero width
    if (in_degree[i] == 0)
    {
      rel_widths[i] = 0.f;
      continue;
    }

    int   outlet = outlet_of[i];
    float max_d = (outlet >= 0) ? max_tree_dist[outlet] : 0.f;

    if (max_d > 1e-6f)
    {
      float norm = 1.f - (dist_to_outlet[i] / max_d);
      rel_widths[i] = std::pow(std::clamp(norm, 0.f, 1.f), exponent);
    }
    else
    {
      rel_widths[i] = 1.f;
    }
  }

  return rel_widths;
}

std::vector<int> Graph::compute_strahler_order() const
{
  const size_t     n = this->size();
  std::vector<int> order(n, 1);

  if (n == 0) return order;

  // build directed downstream graph: edge u -> v if elevation(u) > elevation(v)
  // or tie-broken by node index
  std::vector<std::vector<int>> children(
      n); // upstream contributors to each node
  std::vector<int> out_degree(n, 0);
  std::vector<int> in_degree(n, 0);

  for (const auto &edge : this->edge_list)
  {
    int u = edge.u;
    int v = edge.v;

    if (u < 0 || static_cast<size_t>(u) >= n || v < 0 ||
        static_cast<size_t>(v) >= n || u == v)
      continue;

    float vu = this->points[u].v;
    float vv = this->points[v].v;

    int from = u;
    int to = v;

    if (vu < vv || (vu == vv && u > v))
    {
      from = v;
      to = u;
    }

    // from flows into to
    children[to].push_back(from);
    out_degree[from]++;
    in_degree[to]++;
  }

  // Kahn's algorithm: start with upstream sources (nodes with no incoming
  // upstream flow)
  std::queue<int>  q;
  std::vector<int> in_deg_remaining = in_degree;

  for (size_t i = 0; i < n; ++i)
  {
    if (in_deg_remaining[i] == 0) q.push(static_cast<int>(i));
  }

  std::vector<int> topo_order;
  topo_order.reserve(n);

  // keep map of downstream receiver edges for reverse propagation
  std::vector<std::vector<int>> downstream_neighbors(n);

  for (size_t i = 0; i < n; ++i)
  {
    for (int child : children[i])
      downstream_neighbors[child].push_back(static_cast<int>(i));
  }

  while (!q.empty())
  {
    int u = q.front();
    q.pop();
    topo_order.push_back(u);

    for (int v : downstream_neighbors[u])
    {
      in_deg_remaining[v]--;
      if (in_deg_remaining[v] == 0) q.push(v);
    }
  }

  // if any cycle existed (e.g. flat elevations), append remaining nodes
  if (topo_order.size() < n)
  {
    std::vector<bool> visited(n, false);
    for (int u : topo_order)
      visited[u] = true;

    for (size_t i = 0; i < n; ++i)
      if (!visited[i]) topo_order.push_back(static_cast<int>(i));
  }

  // evaluate Strahler stream order in topological order (upstream to
  // downstream)
  for (int i : topo_order)
  {
    if (children[i].empty())
    {
      order[i] = 1;
      continue;
    }

    int max_child_order = 0;
    int count_max = 0;

    for (int child : children[i])
    {
      int child_order = order[child];

      if (child_order > max_child_order)
      {
        max_child_order = child_order;
        count_max = 1;
      }
      else if (child_order == max_child_order)
      {
        count_max++;
      }
    }

    if (count_max > 1)
      order[i] = max_child_order + 1;
    else
      order[i] = max_child_order;
  }

  return order;
}

void Graph::enforce_downhill_elevation(float min_slope)
{
  const size_t n = this->size();
  if (n == 0 || this->edge_list.empty()) return;

  // determine directed downstream edges from initial elevations
  std::vector<std::vector<int>> downstream_neighbors(n);
  std::vector<int>              in_degree(n, 0);

  for (const auto &edge : this->edge_list)
  {
    int u = edge.u;
    int v = edge.v;

    if (u < 0 || static_cast<size_t>(u) >= n || v < 0 ||
        static_cast<size_t>(v) >= n || u == v)
      continue;

    float vu = this->points[u].v;
    float vv = this->points[v].v;

    int from = u;
    int to = v;

    if (vu < vv || (vu == vv && u > v))
    {
      from = v;
      to = u;
    }

    downstream_neighbors[from].push_back(to);
    in_degree[to]++;
  }

  // Kahn's algorithm: topological sort from upstream headwaters to downstream
  // outlets
  std::queue<int>  q;
  std::vector<int> in_deg_remaining = in_degree;

  for (size_t i = 0; i < n; ++i)
  {
    if (in_deg_remaining[i] == 0) q.push(static_cast<int>(i));
  }

  std::vector<int> topo_order;
  topo_order.reserve(n);

  while (!q.empty())
  {
    int u = q.front();
    q.pop();
    topo_order.push_back(u);

    for (int v : downstream_neighbors[u])
    {
      in_deg_remaining[v]--;
      if (in_deg_remaining[v] == 0) q.push(v);
    }
  }

  // if any cycle existed, append remaining nodes
  if (topo_order.size() < n)
  {
    std::vector<bool> visited(n, false);
    for (int u : topo_order)
      visited[u] = true;

    for (size_t i = 0; i < n; ++i)
      if (!visited[i]) topo_order.push_back(static_cast<int>(i));
  }

  // propagate downstream elevations in topological order
  for (int u : topo_order)
  {
    float zu = this->points[u].v;

    for (int v : downstream_neighbors[u])
    {
      float edge_len = glm::distance(
          glm::vec2(this->points[u].x, this->points[u].y),
          glm::vec2(this->points[v].x, this->points[v].y));
      float dz_min = min_slope * edge_len;

      this->points[v].v = std::min(this->points[v].v, zu - dz_min);
    }
  }
}

size_t Graph::degree(int u) const

{
  if (u < 0 || static_cast<size_t>(u) >= this->adj_list.size()) return 0;
  return this->adj_list[u].size();
}

std::vector<int> Graph::dijkstra(int source_point_index,
                                 int target_point_index) const
{
  if (source_point_index < 0 ||
      source_point_index >= static_cast<int>(this->size()) ||
      target_point_index < 0 ||
      target_point_index >= static_cast<int>(this->size()))
    return {};

  if (source_point_index == target_point_index) return {source_point_index};

  std::vector<float> dist(this->size(), std::numeric_limits<float>::max());
  std::vector<int>   prev(this->size(), -1);

  // --- Dijkstra's algo

  using DistNode = std::pair<float, int>;
  std::priority_queue<DistNode, std::vector<DistNode>, std::greater<DistNode>>
      pq;

  dist[source_point_index] = 0.f;
  pq.push({0.f, source_point_index});

  while (!pq.empty())
  {
    auto [d, u] = pq.top();
    pq.pop();

    if (u == target_point_index) break;
    if (d > dist[u]) continue;

    if (static_cast<size_t>(u) < this->adj_list.size())
    {
      for (const auto &nbr : this->adj_list[u])
      {
        int   k = nbr.target;
        float weight = nbr.weight;
        float alt = dist[u] + weight;
        if (k >= 0 && static_cast<size_t>(k) < this->size() && alt < dist[k])
        {
          dist[k] = alt;
          prev[k] = u;
          pq.push({alt, k});
        }
      }
    }
  }

  if (dist[target_point_index] == std::numeric_limits<float>::max()) return {};

  // --- Backward rebuild the complete path

  std::vector<int> path;
  for (int curr = target_point_index; curr != -1; curr = prev[curr])
  {
    path.push_back(curr);
    if (curr == source_point_index) break;
  }
  std::reverse(path.begin(), path.end());

  return path;
}

bool Graph::empty_edges() const noexcept
{
  return this->edge_list.empty();
}

Edge Graph::get_edge(size_t k) const
{
  if (k < this->edge_list.size()) return this->edge_list[k];
  return {};
}

float Graph::get_edge_length(size_t k) const
{
  if (k >= this->edge_list.size()) return 0.f;
  int u = this->edge_list[k].u;
  int v = this->edge_list[k].v;
  if (u < 0 || u >= static_cast<int>(this->size()) || v < 0 ||
      v >= static_cast<int>(this->size()))
    return 0.f;
  return distance(this->points[u], this->points[v]);
}

std::vector<float> Graph::get_edge_x_pairs() const
{
  std::vector<float> x;
  x.reserve(2 * this->edge_list.size());
  for (const auto &e : this->edge_list)
  {
    if (e.u >= 0 && e.u < static_cast<int>(this->size()) && e.v >= 0 &&
        e.v < static_cast<int>(this->size()))
    {
      x.push_back(this->points[e.u].x);
      x.push_back(this->points[e.v].x);
    }
  }
  return x;
}

std::vector<float> Graph::get_edge_y_pairs() const
{
  std::vector<float> y;
  y.reserve(2 * this->edge_list.size());
  for (const auto &e : this->edge_list)
  {
    if (e.u >= 0 && e.u < static_cast<int>(this->size()) && e.v >= 0 &&
        e.v < static_cast<int>(this->size()))
    {
      y.push_back(this->points[e.u].y);
      y.push_back(this->points[e.v].y);
    }
  }
  return y;
}

const std::vector<Edge> &Graph::get_edges() const noexcept
{
  return this->edge_list;
}

std::vector<float> Graph::get_lengths() const
{
  std::vector<float> lengths;
  lengths.reserve(this->edge_list.size());
  for (size_t k = 0; k < this->edge_list.size(); ++k)
    lengths.push_back(this->get_edge_length(k));
  return lengths;
}

size_t Graph::get_nedges() const noexcept
{
  return this->edge_list.size();
}

Graph Graph::minimum_spanning_tree_prim() const
{
  if (this->empty()) return Graph();

  std::vector<int>   parent(this->size(), -1);
  std::vector<float> key(this->size(), std::numeric_limits<float>::max());
  std::vector<bool>  in_mst(this->size(), false);

  using DistNode = std::pair<float, int>;
  std::priority_queue<DistNode, std::vector<DistNode>, std::greater<DistNode>>
      pq;

  key[0] = 0.f;
  pq.push({0.f, 0});

  while (!pq.empty())
  {
    auto [d, u] = pq.top();
    pq.pop();

    if (in_mst[u] || d > key[u]) continue;
    in_mst[u] = true;

    if (static_cast<size_t>(u) < this->adj_list.size())
    {
      for (const auto &nbr : this->adj_list[u])
      {
        int   v = nbr.target;
        float w = nbr.weight;
        if (v >= 0 && static_cast<size_t>(v) < this->size() && !in_mst[v] &&
            w < key[v])
        {
          key[v] = w;
          parent[v] = u;
          pq.push({w, v});
        }
      }
    }
  }

  // Build output MST graph
  Graph graph(this->points);
  for (size_t i = 1; i < this->size(); ++i)
  {
    if (parent[i] != -1) graph.add_edge(static_cast<int>(i), parent[i], key[i]);
  }

  return graph;
}

const std::vector<Neighbor> &Graph::neighbors(int u) const
{
  if (u < 0 || static_cast<size_t>(u) >= this->adj_list.size())
    return empty_neighbors;
  return this->adj_list[u];
}

size_t Graph::num_edges() const noexcept
{
  return this->edge_list.size();
}

void Graph::print() const
{
  std::cout << "Points:" << std::endl;
  for (size_t k = 0; k < this->size(); k++)
  {
    std::cout << std::setw(6) << k;
    std::cout << std::setw(12) << this->points[k].x;
    std::cout << std::setw(12) << this->points[k].y;
    std::cout << std::setw(12) << this->points[k].v;
    std::cout << std::endl;
  }

  std::cout << "Edges: (index, {pt1, pt2}, weight)" << std::endl;
  for (size_t k = 0; k < this->edge_list.size(); k++)
  {
    std::cout << std::setw(6) << k;
    std::cout << " {" << this->edge_list[k].u << ", ";
    std::cout << this->edge_list[k].v << "} ";
    std::cout << std::setw(12) << this->edge_list[k].weight;
    std::cout << std::endl;
  }
}

Graph Graph::remove_orphan_points() const
{
  Graph            graph_out;
  std::vector<int> new_point_idx(this->size(), -1);

  // Identify connected nodes
  for (size_t u = 0; u < this->size(); ++u)
  {
    if (u < this->adj_list.size() && !this->adj_list[u].empty())
    {
      if (new_point_idx[u] == -1)
      {
        graph_out.push_back(this->points[u]);
        new_point_idx[u] = static_cast<int>(graph_out.size()) - 1;
      }
    }
  }

  // Re-add edges with new indices
  for (const auto &e : this->edge_list)
  {
    if (e.u >= 0 && e.u < static_cast<int>(this->size()) && e.v >= 0 &&
        e.v < static_cast<int>(this->size()))
    {
      int k1 = new_point_idx[e.u];
      int k2 = new_point_idx[e.v];
      if (k1 != -1 && k2 != -1) graph_out.add_edge(k1, k2, e.weight);
    }
  }

  return graph_out;
}

void Graph::set_edge_weight(size_t k, float weight)
{
  if (k >= this->edge_list.size()) return;
  this->edge_list[k].weight = weight;
  int u = this->edge_list[k].u;
  int v = this->edge_list[k].v;
  int edge_idx = static_cast<int>(k);

  if (u >= 0 && static_cast<size_t>(u) < this->adj_list.size())
  {
    for (auto &nbr : this->adj_list[u])
      if (nbr.edge_index == edge_idx) nbr.weight = weight;
  }
  if (v >= 0 && static_cast<size_t>(v) < this->adj_list.size())
  {
    for (auto &nbr : this->adj_list[v])
      if (nbr.edge_index == edge_idx) nbr.weight = weight;
  }
}

void Graph::set_edge_weight(int u, int v, float weight)
{
  if (u < 0 || static_cast<size_t>(u) >= this->adj_list.size()) return;
  for (const auto &nbr : this->adj_list[u])
  {
    if (nbr.target == v)
    {
      this->set_edge_weight(static_cast<size_t>(nbr.edge_index), weight);
      return;
    }
  }
}

void Graph::to_array(Array    &array,
                     glm::vec4 bbox,
                     bool      color_by_edge_weight) const
{
  if (!validate_non_empty(array)) return;

  for (const auto &e : this->edge_list)
  {
    if (e.u >= 0 && e.u < static_cast<int>(this->size()) && e.v >= 0 &&
        e.v < static_cast<int>(this->size()))
    {
      Point p1 = this->points[e.u];
      Point p2 = this->points[e.v];
      if (color_by_edge_weight)
      {
        p1.v = e.weight;
        p2.v = e.weight;
      }
      Path path = Path({p1, p2});
      path.to_array(array, bbox);
    }
  }
}

void Graph::to_array_fractalize(Array        &array,
                                glm::vec4     bbox,
                                int           iterations,
                                std::uint32_t seed,
                                float         sigma,
                                int           orientation,
                                float         persistence) const
{
  if (!validate_non_empty(array)) return;

  // find smallest edge length
  float dmin = std::numeric_limits<float>::max();

  for (size_t k = 0; k < this->edge_list.size(); ++k)
  {
    float dist = this->get_edge_length(k);
    if (dist < dmin && dist > 0.f) dmin = dist;
  }
  if (dmin == std::numeric_limits<float>::max()) dmin = 1.f;

  // fractalize and project to array
  for (const auto &e : this->edge_list)
  {
    if (e.u >= 0 && e.u < static_cast<int>(this->size()) && e.v >= 0 &&
        e.v < static_cast<int>(this->size()))
    {
      Point p1 = this->points[e.u];
      Point p2 = this->points[e.v];
      Path  path = Path({p1, p2});

      path.resample_by_spacing(dmin);
      path = fractalize(path,
                        iterations,
                        seed,
                        sigma,
                        orientation,
                        persistence);
      path.to_array(array, bbox);
    }
  }
}

Array Graph::to_array_sdf(glm::ivec2 shape,
                          glm::vec4  bbox,
                          Array     *p_noise_x,
                          Array     *p_noise_y,
                          glm::vec4  bbox_array) const
{
  (void)bbox;
  return sdf_2d_graph(*this, shape, bbox_array, p_noise_x, p_noise_y, true);
}

void Graph::to_csv(const std::string &fname_xy,
                   const std::string &fname_edges) const
{
  std::ofstream f(fname_xy, std::ios::out);
  for (const auto &p : this->points)
    f << p.x << "," << p.y << "," << p.v << "\n";
  f.close();

  std::ofstream fe(fname_edges, std::ios::out);
  for (const auto &e : this->edge_list)
    fe << e.u << "," << e.v << "," << e.weight << "\n";
  fe.close();
}

void Graph::to_png(const std::string &fname, glm::ivec2 shape) const
{
  if (!validate_shape(shape)) return;

  Array array = Array(shape);
  this->to_array(array, this->get_bbox());
  array.to_png(fname, Cmap::INFERNO, false);
}

} // namespace hmap
