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
#include <limits>
#include <map>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/geometry/graph.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"

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
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  // nodes
  std::vector<float> xp = this->get_edge_x_pairs();
  std::vector<float> yp = this->get_edge_y_pairs();

  for (size_t k = 0; k < xp.size(); k++)
  {
    xp[k] = (xp[k] - bbox.x) / (bbox.y - bbox.x);
    yp[k] = (yp[k] - bbox.z) / (bbox.w - bbox.z);
  }

  // fill heightmap
  auto distance_fct = [&xp, &yp](float x, float y, float)
  {
    float d = std::numeric_limits<float>::max();

    for (size_t i = 0; i < xp.size() - 1; i += 2)
    {
      size_t    j = i + 1;
      glm::vec2 e = {xp[j] - xp[i], yp[j] - yp[i]};
      glm::vec2 w = {x - xp[i], y - yp[i]};
      float     len2 = dot(e, e);
      float     coeff = (len2 > 1e-12f) ? std::clamp(dot(w, e) / len2, 0.f, 1.f)
                                        : 0.f;
      glm::vec2 b = {w.x - e.x * coeff, w.y - e.y * coeff};
      d = std::min(d, dot(b, b));
    }
    return std::sqrt(d);
  };

  Array z = Array(shape);
  fill_array_using_xy_function(z,
                               bbox_array,
                               nullptr,
                               p_noise_x,
                               p_noise_y,
                               nullptr,
                               distance_fct);
  return z;
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
