/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file graph.hpp
 * @brief Definition of the `Graph` class for manipulating 2D graphs.
 *
 * This file contains the definition of the `Graph` class, which extends from
 * the `Cloud` class to handle graph data structures in a 2D space. The `Graph`
 * class supports various functionalities such as graph visualization, edge
 * operations, graph algorithms (e.g., Dijkstra's algorithm, Minimum Spanning
 * Tree), and exporting data to different formats (CSV, PNG).
 *
 * The class includes methods for computing distances, generating visual
 * representations, and manipulating graph edges and nodes. It also provides
 * functionalities for exporting the graph data to CSV and PNG formats, as well
 * as for calculating Signed Distance Functions (SDF) of the graph.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/point.hpp"

namespace hmap
{

class Cloud;

/**
 * @brief Graph edge connecting two node indices with an associated weight.
 */
struct Edge
{
  int   u{0};        ///< First node index
  int   v{0};        ///< Second node index
  float weight{0.f}; ///< Edge weight (cost or length)
};

/**
 * @brief Graph neighbor information for adjacency representation.
 */
struct Neighbor
{
  int   target{0};      ///< Neighbor node index
  float weight{0.f};    ///< Edge weight
  int   edge_index{-1}; ///< Index into the edge list
};

/**
 * @brief Graph class, to manipulate graphs in 2D.
 *
 * This class represents a 2D graph, allowing the creation, manipulation, and
 * analysis of graphs derived from point clouds. It supports operations such as
 * graph construction, traversal, and various geometric analyses. This class
 * inherits from the `Cloud` class, leveraging the functionalities of point
 * clouds while adding graph-specific methods.
 *
 * **Example**
 * @include ex_graph.cpp
 *
 * **Result**
 * @image html ex_graph0.png
 */
class Graph : public Cloud
{
public:
  /**
   * @brief Construct a new Graph object.
   */
  Graph() = default;

  /**
   * @brief Construct a new Graph object based on a cloud of points.
   * @param cloud The cloud of points used to initialize the graph.
   */
  Graph(Cloud cloud);

  /**
   * @brief Construct a new Graph object based on a list of points.
   * @param points The list of points used to initialize the graph.
   */
  Graph(std::vector<Point> points);

  /**
   * @brief Construct a new Graph object based on x and y coordinates.
   * @param x Vector of x coordinates for the points.
   * @param y Vector of y coordinates for the points.
   */
  Graph(std::vector<float> x, std::vector<float> y);

  /**
   * @brief Add an edge connecting two node indices with explicit weight.
   * @param u Node 1 index.
   * @param v Node 2 index.
   * @param weight Edge weight.
   */
  void add_edge(int u, int v, float weight);

  /**
   * @brief Add an edge connecting two node indices with default geometric
   * distance weight.
   * @param u Node 1 index.
   * @param v Node 2 index.
   */
  void add_edge(int u, int v);

  /**
   * @brief Add an edge with explicit weight.
   * @param edge Pair of node indices.
   * @param weight Edge weight.
   */
  void add_edge(glm::ivec2 edge, float weight);

  /**
   * @brief Add an edge with default geometric distance weight.
   * @param edge Pair of node indices.
   */
  void add_edge(glm::ivec2 edge);

  /**
   * @brief Clear all edges from the graph.
   */
  void clear_edges() noexcept;

  /**
   * @brief Get the degree (number of incident edges) of node `u`.
   * @param u Node index.
   * @return Degree of node `u`.
   */
  size_t degree(int u) const;

  /**
   * @brief Return the shortest route between two points using Dijkstra's
   * algorithm.
   *
   * @param  source_point_index Starting point index.
   * @param  target_point_index Ending point index.
   * @return std::vector<int> Path of node indices from source to target.
   *
   * **Example**
   * @include ex_graph_dijkstra.cpp
   */
  std::vector<int> dijkstra(int source_point_index,
                            int target_point_index) const;

  /**
   * @brief Check whether the graph has no edges.
   * @return True if no edges, false otherwise.
   */
  bool empty_edges() const noexcept;

  /**
   * @brief Get an edge by index.
   * @param k Edge index.
   * @return Edge structure.
   */
  Edge get_edge(size_t k) const;

  /**
   * @brief Get the Euclidean length of edge `k`.
   * @param k Edge index.
   * @return float Euclidean length of the edge.
   */
  float get_edge_length(size_t k) const;

  /**
   * @brief Return x coordinates of the edges (as pairs).
   * @return std::vector<float> The x coordinates of the edges.
   */
  std::vector<float> get_edge_x_pairs() const;

  /**
   * @brief Return y coordinates of the edges (as pairs).
   * @return std::vector<float> The y coordinates of the edges.
   */
  std::vector<float> get_edge_y_pairs() const;

  /**
   * @brief Get the edge list.
   * @return Const reference to the edge vector.
   */
  const std::vector<Edge> &get_edges() const noexcept;

  /**
   * @brief Get the Euclidean lengths of all the edges.
   * @return std::vector<float> The lengths of all the edges.
   */
  std::vector<float> get_lengths() const;

  /**
   * @brief Get the number of edges in the graph.
   * @return size_t Number of edges.
   */
  size_t get_nedges() const noexcept;

  /**
   * @brief Generate a Minimum Spanning Tree (MST) using Prim's algorithm.
   * @return Graph The Minimum Spanning Tree (MST) of the graph.
   *
   * **Example**
   * @include ex_graph_minimum_spanning_tree_prim.cpp
   */
  Graph minimum_spanning_tree_prim() const;

  /**
   * @brief Get the neighbors of node `u`.
   * @param u Node index.
   * @return List of Neighbor objects.
   */
  const std::vector<Neighbor> &neighbors(int u) const;

  /**
   * @brief Get the number of edges in the graph.
   * @return size_t Total number of edges.
   */
  size_t num_edges() const noexcept;

  /**
   * @brief Print the graph data to the standard output.
   */
  void print() const;

  /**
   * @brief Remove orphan points (points with no connected edges).
   * @return Graph A new graph object with orphan points removed.
   */
  Graph remove_orphan_points() const;

  /**
   * @brief Set the weight of edge `k`.
   * @param k Edge index.
   * @param weight New weight.
   */
  void set_edge_weight(size_t k, float weight);

  /**
   * @brief Set the weight of an edge between node `u` and node `v`.
   * @param u Node 1 index.
   * @param v Node 2 index.
   * @param weight New weight.
   */
  void set_edge_weight(int u, int v, float weight);

  /**
   * @brief Project the graph to an array and optionally color by edge weight.
   * @param array The input array to project onto.
   * @param bbox Bounding box for projection.
   * @param color_by_edge_weight Color lines by edge weight if true.
   */
  void to_array(Array    &array,
                glm::vec4 bbox,
                bool      color_by_edge_weight = true) const;

  /**
   * @brief Apply fractalization to graph edges and project to an array.
   * @param array Destination array.
   * @param bbox Bounding box.
   * @param iterations Number of fractal iterations.
   * @param seed Random seed.
   * @param sigma Displacement magnitude.
   * @param orientation Displacement orientation.
   * @param persistence Octave persistence.
   */
  void to_array_fractalize(Array        &array,
                           glm::vec4     bbox,
                           int           iterations,
                           std::uint32_t seed,
                           float         sigma = 0.3f,
                           int           orientation = 0,
                           float         persistence = 1.f) const;

  /**
   * @brief Generate an array filled with the Signed Distance Function (SDF) to
   * the graph.
   * @param shape Output array dimensions.
   * @param bbox Bounding box.
   * @param p_noise_x Optional X noise array for domain warping.
   * @param p_noise_y Optional Y noise array for domain warping.
   * @param bbox_array Output array bounding box.
   * @return Array SDF array.
   */
  Array to_array_sdf(glm::ivec2 shape,
                     glm::vec4  bbox,
                     Array     *p_noise_x = nullptr,
                     Array     *p_noise_y = nullptr,
                     glm::vec4  bbox_array = {0.f, 1.f, 0.f, 1.f}) const;

  /**
   * @brief Export graph data to CSV files.
   * @param fname_xy Output CSV for node (x, y, v).
   * @param fname_edges Output CSV for edges (u, v, weight).
   */
  void to_csv(const std::string &fname_xy,
              const std::string &fname_edges) const;

  /**
   * @brief Export the graph as a PNG image file.
   * @param fname Output image path.
   * @param shape Image dimensions in pixels.
   */
  void to_png(const std::string &fname, glm::ivec2 shape = {512, 512}) const;

private:
  std::vector<Edge>                  edge_list;
  std::vector<std::vector<Neighbor>> adj_list;

  static const std::vector<Neighbor> empty_neighbors;

  void ensure_adj_capacity(size_t node_count);
};

} // namespace hmap