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

TEST(GraphTest, DistanceToOutlet)
{
  // Y-junction:
  // Node 0 (headwater 1, elev 1.0) -> Node 2 (junction, elev 0.5)
  // Node 1 (headwater 2, elev 1.0) -> Node 2 (junction, elev 0.5)
  // Node 2 (junction, elev 0.5)    -> Node 3 (outlet, elev 0.0)
  Graph g;
  g.push_back(Point(0.5f, 0.9f, 1.0f)); // 0
  g.push_back(Point(0.9f, 0.5f, 1.0f)); // 1
  g.push_back(Point(0.5f, 0.5f, 0.5f)); // 2
  g.push_back(Point(0.5f, 0.1f, 0.0f)); // 3

  g.add_edge(0, 2);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  std::vector<float> dist = g.compute_distance_to_outlet();
  ASSERT_EQ(dist.size(), 4u);

  // outlet has distance 0
  EXPECT_FLOAT_EQ(dist[3], 0.0f);

  // node 2 is 0.4 from outlet
  EXPECT_NEAR(dist[2], 0.4f, 1e-5f);

  // node 0 is 0.4 from node 2, so 0.8 from outlet
  EXPECT_NEAR(dist[0], 0.8f, 1e-5f);

  // node 1 is 0.4 from node 2, so 0.8 from outlet
  EXPECT_NEAR(dist[1], 0.8f, 1e-5f);
}

TEST(GraphTest, StrahlerOrder)
{
  // Y-junction:
  // Node 0 (headwater 1, elev 1.0) -> Node 2 (junction, elev 0.5)
  // Node 1 (headwater 2, elev 1.0) -> Node 2 (junction, elev 0.5)
  // Node 2 (junction, elev 0.5)    -> Node 3 (outlet, elev 0.0)
  Graph g;
  g.push_back(Point(0.0f, 1.0f, 1.0f)); // 0
  g.push_back(Point(1.0f, 1.0f, 1.0f)); // 1
  g.push_back(Point(0.5f, 0.5f, 0.5f)); // 2
  g.push_back(Point(0.5f, 0.0f, 0.0f)); // 3

  g.add_edge(0, 2);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  std::vector<int> strahler = g.compute_strahler_order();
  ASSERT_EQ(strahler.size(), 4u);

  EXPECT_EQ(strahler[0], 1);
  EXPECT_EQ(strahler[1], 1);
  EXPECT_EQ(strahler[2], 2); // 1 + 1 confluent -> order 2
  EXPECT_EQ(strahler[3], 2); // single order 2 branch entering -> order 2
}

TEST(GraphTest, EnforceDownhillElevation)
{
  // Test with flat spots / uphill anomalies:
  // Node 0 (headwater 1, elev 1.0) -> Node 2 (junction, elev 0.5) -> Node 3
  // (outlet, elev 0.52 - slight uphill error) Node 1 (headwater 2, elev 0.8) ->
  // Node 2 (junction, elev 0.5)
  Graph g;
  g.push_back(Point(0.0f, 1.0f, 1.0f)); // 0: (0, 1) elev 1.0
  g.push_back(Point(1.0f, 1.0f, 0.8f)); // 1: (1, 1) elev 0.8
  g.push_back(Point(0.5f, 0.5f, 0.5f)); // 2: (0.5, 0.5) elev 0.5
  g.push_back(Point(
      0.5f,
      0.0f,
      0.52f)); // 3: (0.5, 0.0) elev 0.52 (small anomaly relative to junction)

  // In graph construction from stream paths, edges are directed from path order
  // or downstream flow Here 0 -> 2, 1 -> 2, and 2 -> 3
  g.add_edge(0, 2);
  g.add_edge(1, 2);
  // Edge 2-3 has u=2, v=3. Since v2=0.5 and v3=0.52, if based solely on initial
  // v, 3 would flow to 2! But with 0.48 for node 3 and a flat spot v=0.5 at
  // junction: Suppose node 2 is flat with node 0 (elev 1.0):
  Graph g2;
  g2.push_back(Point(0.0f, 1.0f, 1.0f)); // 0: elev 1.0
  g2.push_back(Point(1.0f, 1.0f, 1.0f)); // 1: elev 1.0
  g2.push_back(Point(0.5f, 0.5f, 1.0f)); // 2: flat spot (elev 1.0)
  g2.push_back(Point(0.5f, 0.0f, 0.2f)); // 3: outlet (elev 0.2)

  g2.add_edge(0, 2);
  g2.add_edge(1, 2);
  g2.add_edge(2, 3);

  float min_slope = 0.1f;
  g2.enforce_downhill_elevation(min_slope);

  // node 0 initial elevation is 1.0
  EXPECT_FLOAT_EQ(g2[0].v, 1.0f);
  EXPECT_FLOAT_EQ(g2[1].v, 1.0f);

  // node 2 should be decreased below 1.0
  float dist_02 = std::sqrt(0.5f * 0.5f + 0.5f * 0.5f);
  EXPECT_LE(g2[2].v, 1.0f - min_slope * dist_02);
  EXPECT_LT(g2[2].v, 1.0f);

  // node 3 should be strictly lower than node 2
  float dist_23 = 0.5f;
  EXPECT_LE(g2[3].v, g2[2].v - min_slope * dist_23);
  EXPECT_LT(g2[3].v, g2[2].v);
}

TEST(GraphTest, FlowAccumulation)
{
  // A Y-shaped confluence network:
  // Node 0 (0, 1) and Node 1 (1, 1) meet at Node 2 (0.5, 0.5), which flows to
  // Node 3 (0.5, 0.0)
  Graph g;
  g.push_back(Point(0.0f, 1.0f, 1.0f)); // 0: headwater 1
  g.push_back(Point(1.0f, 1.0f, 0.9f)); // 1: headwater 2
  g.push_back(Point(0.5f, 0.5f, 0.5f)); // 2: confluence
  g.push_back(Point(0.5f, 0.0f, 0.1f)); // 3: outlet

  g.add_edge(0, 2);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  std::vector<float> accum = g.compute_flow_accumulation();
  ASSERT_EQ(accum.size(), 4u);

  // headwaters have 0 accumulation
  EXPECT_FLOAT_EQ(accum[0], 0.f);
  EXPECT_FLOAT_EQ(accum[1], 0.f);

  // node 2 accumulates length from edge (0->2) and edge (1->2)
  float len_02 = std::sqrt(0.5f * 0.5f + 0.5f * 0.5f);
  float len_12 = std::sqrt(0.5f * 0.5f + 0.5f * 0.5f);
  EXPECT_NEAR(accum[2], len_02 + len_12, 1e-4f);

  // node 3 accumulates accum[2] plus edge (2->3)
  float len_23 = 0.5f;
  EXPECT_NEAR(accum[3], accum[2] + len_23, 1e-4f);
  EXPECT_GT(accum[3], accum[2]);
}
