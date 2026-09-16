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

std::vector<int> Graph::dijkstra(int source_point_index, int target_point_index)
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

    for (int k : this->connectivity[u])
    {
      auto it = this->adjacency_matrix.find({u, k});
      if (it == this->adjacency_matrix.end()) continue;

      float weight = it->second;
      float alt = dist[u] + weight;
      if (alt < dist[k])
      {
        dist[k] = alt;
        prev[k] = u;
        pq.push({alt, k});
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

void Graph::add_edge(std::vector<int> edge, float weight)
{
  this->edges.push_back(edge);
  this->weights.push_back(weight);
}

void Graph::add_edge(std::vector<int> edge)
{
  this->edges.push_back(edge);
  this->weights.push_back(this->get_edge_length((int)this->get_nedges() - 1));
}

float Graph::get_edge_length(int k)
{
  return distance(this->points[this->edges[k][0]],
                  this->points[this->edges[k][1]]);
}

std::vector<float> Graph::get_lengths()
{
  std::vector<float> lengths = {};
  for (size_t k = 0; k < this->get_nedges(); k++)
  {
    lengths.push_back(distance(this->points[this->edges[k][0]],
                               this->points[this->edges[k][1]]));
  }
  return lengths;
}

std::vector<float> Graph::get_edge_x_pairs()
{
  std::vector<float> x;
  x.reserve(2 * this->get_nedges());
  for (auto &e : this->edges)
  {
    x.push_back(this->points[e[0]].x);
    x.push_back(this->points[e[1]].x);
  }
  return x;
}

std::vector<float> Graph::get_edge_y_pairs()
{
  std::vector<float> y;
  y.reserve(2 * this->get_nedges());
  for (auto &e : this->edges)
  {
    y.push_back(this->points[e[0]].y);
    y.push_back(this->points[e[1]].y);
  }
  return y;
}

size_t Graph::get_nedges()
{
  return this->edges.size();
}

Graph Graph::minimum_spanning_tree_prim()
{
  std::vector<int>   parent(this->size());
  std::vector<float> key(this->size());
  std::vector<bool>  is_point_in_mst(this->size());

  for (size_t i = 0; i < this->size(); i++)
  {
    key[i] = std::numeric_limits<float>::max();
    is_point_in_mst[i] = false;
  }

  // starting point
  key[0] = 0.f;
  parent[0] = -1;

  for (size_t i = 0; i < this->size() - 1; i++)
  {
    // find point with smallest 'key' while not being in the MS tree
    int   k = 0;
    float key_max = std::numeric_limits<float>::max();
    for (size_t p = 0; p < key.size(); p++)
      if ((key[p] < key_max) and (is_point_in_mst[p] == false))
      {
        key_max = key[p];
        k = (int)p;
      }

    is_point_in_mst[k] = true;

    for (size_t p = 0; p < this->size(); p++)
    {
      if ((this->adjacency_matrix[{k, p}] > 0.f) and
          (is_point_in_mst[p] == false) and
          (this->adjacency_matrix[{k, p}] < key[p]))
      {
        parent[p] = k;
        key[p] = this->adjacency_matrix[{k, p}];
      }
    }
  }

  // build output graph
  Graph graph = Graph(this->points);
  for (size_t i = 1; i < this->size(); i++)
    graph.add_edge({(int)i, parent[i]});

  return graph;
}

void Graph::print()
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
  for (size_t k = 0; k < this->get_nedges(); k++)
  {
    std::cout << std::setw(6) << k;
    std::cout << " {" << this->edges[k][0] << ", ";
    std::cout << this->edges[k][1] << "} ";
    std::cout << std::setw(12) << this->weights[k];
    std::cout << std::endl;
  }
}

Graph Graph::remove_orphan_points()
{
  Graph            graph_out = Graph();
  std::vector<int> new_point_idx(this->size());

  // fill vector with '-1' to keep track of which points have already
  // been added
  std::fill(new_point_idx.begin(), new_point_idx.end(), -1);

  this->update_connectivity();

  for (size_t k = 0; k < this->size(); k++)
  {
    if (this->connectivity[k].size() > 0)
    {
      // current point is connected to at least one other node =>
      // add it
      if (new_point_idx[k] == -1)
      {
        graph_out.add_point(this->points[k]);
        new_point_idx[k] = (int)graph_out.size() - 1;
      }

      for (size_t r = 0; r < this->connectivity[k].size(); r++)
      {
        int j = this->connectivity[k][r];
        if ((j > (int)k) and (new_point_idx[j] == -1))
        {
          graph_out.add_point(this->points[j]);
          new_point_idx[j] = (int)graph_out.size() - 1;
        }
      }
    }
  }

  // rebuild connectivity
  for (size_t k = 0; k < this->get_nedges(); k++)
  {
    int k1 = new_point_idx[this->edges[k][0]];
    int k2 = new_point_idx[this->edges[k][1]];
    graph_out.add_edge({k1, k2}, this->weights[k]);
  }

  return graph_out;
}

void Graph::to_array(Array &array, glm::vec4 bbox, bool color_by_edge_weight)
{
  if (!validate_non_empty(array)) return;

  if (color_by_edge_weight)
    for (std::size_t k = 0; k < this->get_nedges(); k++)
    {
      Point p1 = this->points[this->edges[k][0]];
      Point p2 = this->points[this->edges[k][1]];
      p1.v = this->weights[this->edges[k][0]];
      p2.v = this->weights[this->edges[k][1]];
      Path path = Path({p1, p2});
      path.to_array(array, bbox);
    }
  else
    for (std::size_t k = 0; k < this->get_nedges(); k++)
    {
      Point p1 = this->points[this->edges[k][0]];
      Point p2 = this->points[this->edges[k][1]];
      Path  path = Path({p1, p2});
      path.to_array(array, bbox);
    }
}

void Graph::to_array_fractalize(Array        &array,
                                glm::vec4     bbox,
                                int           iterations,
                                std::uint32_t seed,
                                float         sigma,
                                int           orientation,
                                float         persistence)
{
  if (!validate_non_empty(array)) return;

  // find smallest edge length
  float dmin = std::numeric_limits<float>::max();

  for (size_t k = 0; k < this->get_nedges(); k++)
  {
    float dist = this->get_edge_length((int)k);
    if (dist < dmin) dmin = dist;
  }

  // fractalize and project to array
  for (std::size_t k = 0; k < this->get_nedges(); k++)
  {
    Point p1 = this->points[this->edges[k][0]];
    Point p2 = this->points[this->edges[k][1]];
    Path  path = Path({p1, p2});

    path.resample_by_spacing(dmin);
    path = fractalize(path, iterations, seed, sigma, orientation, persistence);
    path.to_array(array, bbox);
  }
}

Array Graph::to_array_sdf(glm::ivec2 shape,
                          glm::vec4  bbox,
                          Array     *p_noise_x,
                          Array     *p_noise_y,
                          glm::vec4  bbox_array)
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
      float     coeff = std::clamp(dot(w, e) / dot(e, e), 0.f, 1.f);
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

void Graph::to_csv(std::string fname_xy, std::string fname_adjacency)
{
  std::fstream f;

  f.open(fname_xy, std::ios::out);
  for (auto &p : this->points)
    f << p.x << "," << p.y << "," << p.v << std::endl;
  f.close();

  f.open(fname_adjacency, std::ios::out);
  for (int i = 0; i < (int)this->size(); i++)
  {
    for (int j = 0; j < (int)this->size(); j++)
    {
      float v = 0.f;
      if (this->adjacency_matrix.count({i, j}))
        v = this->adjacency_matrix[{i, j}];
      f << v;

      if (j < (int)this->size() - 1) f << ",";
    }
    f << std::endl;
  }

  f.close();
}

void Graph::to_png(std::string fname, glm::ivec2 shape)
{
  if (!validate_shape(shape)) return;

  Array array = Array(shape);
  this->to_array(array, this->get_bbox());
  array.to_png(fname, Cmap::INFERNO, false);
}

void Graph::update_adjacency_matrix()
{
  this->adjacency_matrix.clear();

  // fill matrix
  for (std::size_t k = 0; k < this->get_nedges(); k++)
  {
    this->adjacency_matrix[{this->edges[k][0], this->edges[k][1]}] =
        this->weights[k];
    this->adjacency_matrix[{this->edges[k][1], this->edges[k][0]}] =
        this->weights[k];
  }
}

void Graph::update_connectivity()
{
  std::vector<std::vector<int>> nbrs(this->size());

  for (std::size_t k = 0; k < this->get_nedges(); k++)
  {
    nbrs[this->edges[k][0]].push_back(this->edges[k][1]);
    nbrs[this->edges[k][1]].push_back(this->edges[k][0]);
  }

  this->connectivity = nbrs;
}

} // namespace hmap
