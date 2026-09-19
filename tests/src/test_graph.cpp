#include <cmath>
#include <vector>

#include "highmap/geometry/graph.hpp"
#include "highmap/geometry/point.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(GraphTest, ConstructionAndEdgeManagement)
{
  Graph g;
  EXPECT_TRUE(g.empty());
  EXPECT_TRUE(g.empty_edges());
  EXPECT_EQ(g.num_edges(), 0u);
  EXPECT_EQ(g.get_nedges(), 0u);

  g.push_back(Point(0.f, 0.f));
  g.push_back(Point(1.f, 0.f));
  g.push_back(Point(1.f, 1.f));
  EXPECT_EQ(g.size(), 3u);

  g.add_edge(0, 1);
  EXPECT_EQ(g.num_edges(), 1u);
  EXPECT_FLOAT_EQ(g.get_edge(0).weight, 1.f);
  EXPECT_FLOAT_EQ(g.get_edge_length(0), 1.f);

  g.add_edge(glm::ivec2(1, 2), 2.5f);
  EXPECT_EQ(g.num_edges(), 2u);
  EXPECT_FLOAT_EQ(g.get_edge(1).weight, 2.5f);
  EXPECT_FLOAT_EQ(g.get_edge_length(1), 1.f);

  EXPECT_EQ(g.degree(0), 1u);
  EXPECT_EQ(g.degree(1), 2u);
  EXPECT_EQ(g.degree(2), 1u);

  const auto &nbrs1 = g.neighbors(1);
  ASSERT_EQ(nbrs1.size(), 2u);
  EXPECT_EQ(nbrs1[0].target, 0);
  EXPECT_EQ(nbrs1[1].target, 2);

  g.set_edge_weight(0, 10.f);
  EXPECT_FLOAT_EQ(g.get_edge(0).weight, 10.f);
  EXPECT_FLOAT_EQ(g.neighbors(1)[0].weight, 10.f);

  g.clear_edges();
  EXPECT_TRUE(g.empty_edges());
  EXPECT_EQ(g.degree(1), 0u);
}

TEST(GraphTest, DijkstraShortestPath)
{
  // 0 - (1) - 1 - (2) - 2
  // |                   |
  // +------- (10) ------+
  Graph g;
  g.push_back(Point(0.f, 0.f));
  g.push_back(Point(1.f, 0.f));
  g.push_back(Point(2.f, 0.f));

  g.add_edge(0, 1, 1.f);
  g.add_edge(1, 2, 2.f);
  g.add_edge(0, 2, 10.f);

  std::vector<int> path = g.dijkstra(0, 2);
  ASSERT_EQ(path.size(), 3u);
  EXPECT_EQ(path[0], 0);
  EXPECT_EQ(path[1], 1);
  EXPECT_EQ(path[2], 2);

  // Shortest path with updated direct edge weight
  g.set_edge_weight(0, 2, 0.5f);
  path = g.dijkstra(0, 2);
  ASSERT_EQ(path.size(), 2u);
  EXPECT_EQ(path[0], 0);
  EXPECT_EQ(path[1], 2);
}

TEST(GraphTest, MinimumSpanningTreePrim)
{
  // Triangle graph: (0, 1)=1, (1, 2)=2, (0, 2)=3
  Graph g;
  g.push_back(Point(0.f, 0.f));
  g.push_back(Point(1.f, 0.f));
  g.push_back(Point(0.f, 1.f));

  g.add_edge(0, 1, 1.f);
  g.add_edge(1, 2, 2.f);
  g.add_edge(0, 2, 3.f);

  Graph mst = g.minimum_spanning_tree_prim();
  EXPECT_EQ(mst.size(), 3u);
  EXPECT_EQ(mst.num_edges(), 2u);

  float total_weight = 0.f;
  for (const auto &e : mst.get_edges())
    total_weight += e.weight;
  EXPECT_FLOAT_EQ(total_weight, 3.f); // 1 + 2
}

TEST(GraphTest, RemoveOrphanPoints)
{
  Graph g;
  g.push_back(Point(0.f, 0.f)); // node 0 (orphan)
  g.push_back(Point(1.f, 1.f)); // node 1
  g.push_back(Point(2.f, 2.f)); // node 2
  g.push_back(Point(3.f, 3.f)); // node 3 (orphan)

  g.add_edge(1, 2, 5.f);

  Graph cleaned = g.remove_orphan_points();
  EXPECT_EQ(cleaned.size(), 2u);
  EXPECT_EQ(cleaned.num_edges(), 1u);
  EXPECT_FLOAT_EQ(cleaned[0].x, 1.f);
  EXPECT_FLOAT_EQ(cleaned[0].y, 1.f);
  EXPECT_FLOAT_EQ(cleaned[1].x, 2.f);
  EXPECT_FLOAT_EQ(cleaned[1].y, 2.f);
  EXPECT_EQ(cleaned.get_edge(0).u, 0);
  EXPECT_EQ(cleaned.get_edge(0).v, 1);
  EXPECT_FLOAT_EQ(cleaned.get_edge(0).weight, 5.f);
}

TEST(GraphTest, EdgeCoordinatesAndLengths)
{
  Graph g;
  g.push_back(Point(0.f, 0.f));
  g.push_back(Point(3.f, 4.f));

  g.add_edge(0, 1);

  auto x_pairs = g.get_edge_x_pairs();
  auto y_pairs = g.get_edge_y_pairs();
  auto lengths = g.get_lengths();

  ASSERT_EQ(x_pairs.size(), 2u);
  ASSERT_EQ(y_pairs.size(), 2u);
  ASSERT_EQ(lengths.size(), 1u);

  EXPECT_FLOAT_EQ(x_pairs[0], 0.f);
  EXPECT_FLOAT_EQ(x_pairs[1], 3.f);
  EXPECT_FLOAT_EQ(y_pairs[0], 0.f);
  EXPECT_FLOAT_EQ(y_pairs[1], 4.f);
  EXPECT_FLOAT_EQ(lengths[0], 5.f);
}
