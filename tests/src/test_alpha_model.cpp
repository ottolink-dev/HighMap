#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/geometry/graph.hpp"
#include "highmap/roads.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(AlphaModelTest, InvalidInputs)
{
  Array     z = Array(glm::ivec2(16, 16), 0.f);
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};

  // empty cities
  Graph g1 = generate_network_alpha_model({}, {}, {}, bbox, z, 0);
  EXPECT_EQ(g1.size(), 0);

  // single city (< 2)
  Graph g2 = generate_network_alpha_model({0.5f}, {0.5f}, {1.f}, bbox, z, 0);
  EXPECT_EQ(g2.size(), 0);

  // mismatched vector sizes
  Graph g3 = generate_network_alpha_model({0.2f, 0.8f},
                                          {0.2f},
                                          {1.f, 1.f},
                                          bbox,
                                          z,
                                          0);
  EXPECT_EQ(g3.size(), 0);

  // empty array
  Array z_empty;
  Graph g4 = generate_network_alpha_model({0.2f, 0.8f},
                                          {0.2f, 0.8f},
                                          {1.f, 1.f},
                                          bbox,
                                          z_empty,
                                          0);
  EXPECT_EQ(g4.size(), 0);
}

TEST(AlphaModelTest, BasicNetworkGeneration)
{
  glm::ivec2 shape = {64, 64};
  Array      z = Array(shape, 0.f);
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  std::vector<float> xc = {0.2f, 0.8f, 0.5f};
  std::vector<float> yc = {0.2f, 0.8f, 0.5f};
  std::vector<float> size = {10.f, 10.f, 5.f};

  Graph net =
      generate_network_alpha_model(xc, yc, size, bbox, z, 42, 0.7f, 100);
  EXPECT_GT(net.size(), 0);
  EXPECT_GT(net.get_nedges(), 0);
}

TEST(GraphDijkstraTest, PathWithIntermediateZero)
{
  // Build a linear graph: 1 -> 0 -> 2
  Graph g;
  g.push_back(Point(0.f, 0.f, 0.f));  // node 0
  g.push_back(Point(-1.f, 0.f, 0.f)); // node 1
  g.push_back(Point(1.f, 0.f, 0.f));  // node 2

  g.add_edge({1, 0}, 1.f);
  g.add_edge({0, 2}, 1.f);
  g.update_adjacency_matrix();
  g.update_connectivity();

  std::vector<int> path = g.dijkstra(1, 2);
  ASSERT_EQ(path.size(), 3);
  EXPECT_EQ(path[0], 1);
  EXPECT_EQ(path[1], 0);
  EXPECT_EQ(path[2], 2);
}

TEST(GraphDijkstraTest, SameSourceAndTarget)
{
  Graph g;
  g.push_back(Point(0.f, 0.f, 0.f));
  g.push_back(Point(1.f, 0.f, 0.f));
  g.add_edge({0, 1}, 1.f);
  g.update_adjacency_matrix();
  g.update_connectivity();

  std::vector<int> path = g.dijkstra(0, 0);
  ASSERT_EQ(path.size(), 1);
  EXPECT_EQ(path[0], 0);
}

TEST(GraphDijkstraTest, DisconnectedGraph)
{
  Graph g;
  g.push_back(Point(0.f, 0.f, 0.f)); // node 0
  g.push_back(Point(1.f, 0.f, 0.f)); // node 1
  g.update_adjacency_matrix();
  g.update_connectivity();

  std::vector<int> path = g.dijkstra(0, 1);
  EXPECT_TRUE(path.empty());
}
