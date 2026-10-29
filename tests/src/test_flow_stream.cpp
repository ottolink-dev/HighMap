#include <cmath>
#include <vector>

#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(FlowStreamMesh, StreamToBoundaryVertex)
{
  // inclined plane from left (high) to right (low)
  std::vector<glm::vec3> pts = {{0.0f, 0.0f, 1.0f},
                                {0.0f, 1.0f, 1.0f},
                                {1.0f, 0.0f, 0.0f},
                                {1.0f, 1.0f, 0.0f},
                                {0.5f, 0.5f, 0.5f}};

  hmap::TerrainTriMesh mesh(pts);

  // start from top-left (vertex 1)
  hmap::Path path = hmap::flow_stream(mesh, size_t(1));

  EXPECT_GE(path.size(), 2u);

  // path must end at a boundary outlet
  auto last_pt = path[path.size() - 1];
  EXPECT_NEAR(last_pt.x, 1.0f, 1e-4f);
  EXPECT_NEAR(last_pt.v, 0.0f, 1e-4f);
}

TEST(FlowStreamMesh, StreamFromCoordinates)
{
  std::vector<glm::vec3> pts = {{0.0f, 0.0f, 1.0f},
                                {0.0f, 1.0f, 1.0f},
                                {1.0f, 0.0f, 0.0f},
                                {1.0f, 1.0f, 0.0f},
                                {0.5f, 0.5f, 0.5f}};

  hmap::TerrainTriMesh mesh(pts);

  // start near (0.05, 0.95)
  hmap::Path path = hmap::flow_stream(mesh, glm::vec2(0.05f, 0.95f));

  EXPECT_GE(path.size(), 2u);

  auto last_pt = path[path.size() - 1];
  EXPECT_NEAR(last_pt.x, 1.0f, 1e-4f);
}

TEST(FlowStreamMesh, StreamFromHeightmapMesh)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);
  hmap::remap(z);

  hmap::TerrainTriMesh mesh =
      hmap::generate_terrain_tri_mesh_from_heightmap_random(z, 256, 42);

  hmap::Path path = hmap::flow_stream(mesh, glm::vec2(0.5f, 0.5f));

  EXPECT_GE(path.size(), 2u);

  for (size_t i = 0; i < path.size(); ++i)
  {
    const auto &pt = path[i];
    EXPECT_GE(pt.x, 0.f);
    EXPECT_LE(pt.x, 1.f);
    EXPECT_GE(pt.y, 0.f);
    EXPECT_LE(pt.y, 1.f);
  }
}

TEST(FlowStreamVirtualArray, StreamOnVirtualArray)
{
  glm::ivec2        shape = {64, 64};
  glm::ivec2        tile_shape = {32, 32};
  int               halo = 8;
  hmap::ComputeMode cm;

  hmap::Array z_arr = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {2.f, 2.f},
                                      42);
  hmap::remap(z_arr);

  hmap::VirtualArray z_va(shape, tile_shape, halo, hmap::StorageMode::VA_RAM);
  z_va.from_array(z_arr, cm);

  hmap::Path path = hmap::va::flow_stream(cm,
                                          z_va,
                                          glm::vec2(0.4f, 0.6f),
                                          42,
                                          256);

  EXPECT_GE(path.size(), 2u);

  for (size_t i = 0; i < path.size(); ++i)
  {
    const auto &pt = path[i];
    EXPECT_GE(pt.x, 0.f);
    EXPECT_LE(pt.x, 1.f);
    EXPECT_GE(pt.y, 0.f);
    EXPECT_LE(pt.y, 1.f);
  }
}

TEST(FlowStreamArray, FavorExistingStreams)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  123);
  hmap::remap(z);

  // generate first stream
  glm::ivec2 ij_start1(50, 64);
  hmap::Path stream1 = hmap::flow_stream(z, ij_start1);

  ASSERT_GE(stream1.size(), 2u);

  // start stream2 near stream1 (offset by 1 grid cell)
  size_t     mid_idx = stream1.size() / 2;
  glm::ivec2 ij_start2(
      std::clamp(
          static_cast<int>(std::round(stream1[mid_idx].x * (shape.x - 1.f))) +
              1,
          0,
          shape.x - 1),
      std::clamp(
          static_cast<int>(std::round(stream1[mid_idx].y * (shape.y - 1.f))) +
              1,
          0,
          shape.y - 1));

  hmap::Path stream2_independent = hmap::flow_stream(z, ij_start2);
  ASSERT_GE(stream2_independent.size(), 2u);

  // generate second stream with favoring existing streams
  hmap::Path stream2_favored = hmap::flow_stream(z,
                                                 ij_start2,
                                                 0.5f,
                                                 2.f,
                                                 100.f,
                                                 {stream1});

  ASSERT_GE(stream2_favored.size(), 2u);

  // end point of favored stream should match the exit point of stream1
  const auto &end1 = stream1[stream1.size() - 1];
  const auto &end_favored = stream2_favored[stream2_favored.size() - 1];

  EXPECT_NEAR(end1.x, end_favored.x, 1e-4f);
  EXPECT_NEAR(end1.y, end_favored.y, 1e-4f);
}

TEST(FlowStreamMesh, FavorExistingStreams)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);
  hmap::remap(z);

  hmap::TerrainTriMesh mesh =
      hmap::generate_terrain_tri_mesh_from_heightmap_random(z, 512, 42);

  glm::vec2  start1(0.4f, 0.5f);
  hmap::Path stream1 = hmap::flow_stream(mesh, start1);

  ASSERT_GE(stream1.size(), 2u);

  glm::vec2  start2(0.45f, 0.55f);
  hmap::Path stream2_favored =
      hmap::flow_stream(mesh, start2, 0.5f, 2.f, 100.f, false, {stream1});

  ASSERT_GE(stream2_favored.size(), 2u);

  const auto &end1 = stream1[stream1.size() - 1];
  const auto &end_favored = stream2_favored[stream2_favored.size() - 1];

  EXPECT_NEAR(end1.x, end_favored.x, 1e-3f);
  EXPECT_NEAR(end1.y, end_favored.y, 1e-3f);
}

TEST(FlowStreamMesh, NetworkReachBoundaryAndDecreasingElevation)
{
  glm::ivec2 shape = {128, 128};
  int        seed = 42;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {3.f, 3.f},
                                  seed);
  z += hmap::circus(shape);
  hmap::remap(z);

  hmap::TerrainTriMesh mesh =
      hmap::generate_terrain_tri_mesh_from_heightmap_random(z, 2048, seed);

  std::vector<glm::vec2> start_points = {
      {0.35f, 0.40f},
      {0.30f, 0.60f},
      {0.65f, 0.65f},
      {0.45f, 0.70f},
  };

  // 1. Individual stream reachability to boundary
  for (const auto &start : start_points)
  {
    hmap::Path path = hmap::flow_stream(mesh, start);
    ASSERT_GE(path.size(), 2u);

    const auto &last_pt = path.back();
    bool        on_boundary = (last_pt.x <= 1e-3f || last_pt.x >= 1.f - 1e-3f ||
                        last_pt.y <= 1e-3f || last_pt.y >= 1.f - 1e-3f);
    EXPECT_TRUE(on_boundary);
  }

  // 2. Sequential network generation with favoring and intersection truncation
  std::vector<hmap::Path> streams;

  for (const auto &start : start_points)
  {
    hmap::Path stream = hmap::flow_stream(mesh,
                                          start,
                                          0.95f,
                                          2.f,
                                          100.f,
                                          true,
                                          streams,
                                          0.03f,
                                          0.4f,
                                          true);

    if (!stream.empty()) streams.push_back(stream);
  }

  ASSERT_EQ(streams.size(), start_points.size());

  // 3. Convert to graph and verify network connectivity
  hmap::Graph stream_graph = hmap::hydrology_paths_to_graph(streams, 1e-3f);
  ASSERT_GT(stream_graph.size(), 0u);
  ASSERT_GT(stream_graph.num_edges(), 0u);

  // Enforce downhill elevation
  float min_slope = 0.005f;
  stream_graph.enforce_downhill_elevation(min_slope);

  // Build adjacency
  std::vector<std::vector<int>> adj(stream_graph.size());

  for (size_t e = 0; e < stream_graph.num_edges(); ++e)
  {
    auto edge = stream_graph.get_edge(e);
    adj[edge.u].push_back(edge.v);
    adj[edge.v].push_back(edge.u);
  }

  // Verify that tracing downstream from every stream starting point
  // strictly decreases in elevation and terminates at a boundary outlet node
  for (size_t s = 0; s < streams.size(); ++s)
  {
    int start_node = -1;

    for (size_t i = 0; i < stream_graph.size(); ++i)
    {
      if (std::abs(stream_graph[i].x - streams[s][0].x) < 1e-4f &&
          std::abs(stream_graph[i].y - streams[s][0].y) < 1e-4f)
      {
        start_node = static_cast<int>(i);
        break;
      }
    }

    ASSERT_NE(start_node, -1);

    int  curr = start_node;
    bool reached_boundary = false;
    int  steps = 0;

    while (steps < 300)
    {
      if (stream_graph[curr].x <= 1e-3f ||
          stream_graph[curr].x >= 1.f - 1e-3f ||
          stream_graph[curr].y <= 1e-3f || stream_graph[curr].y >= 1.f - 1e-3f)
      {
        reached_boundary = true;
        break;
      }

      int   next_node = -1;
      float best_z = stream_graph[curr].v;

      for (int nbr : adj[curr])
      {
        if (stream_graph[nbr].v < best_z)
        {
          best_z = stream_graph[nbr].v;
          next_node = nbr;
        }
      }

      if (next_node == -1) break;

      // strictly decreasing elevation check
      EXPECT_LT(stream_graph[next_node].v, stream_graph[curr].v);

      curr = next_node;
      steps++;
    }

    EXPECT_TRUE(reached_boundary);
  }

  // Verify compute_distance_to_outlet on the stream network
  std::vector<float> dist_to_outlet = stream_graph.compute_distance_to_outlet();

  for (size_t s = 0; s < streams.size(); ++s)
  {
    for (size_t i = 0; i < stream_graph.size(); ++i)
    {
      if (std::abs(stream_graph[i].x - streams[s][0].x) < 1e-4f &&
          std::abs(stream_graph[i].y - streams[s][0].y) < 1e-4f)
      {
        EXPECT_GT(dist_to_outlet[i], 0.f);
        break;
      }
    }
  }
}

TEST(FlowStream, TruncateAtIntersection)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);
  hmap::remap(z);

  // straight horizontal line across the domain
  hmap::Path existing(
      {hmap::Point(0.1f, 0.5f, 0.5f), hmap::Point(0.9f, 0.5f, 0.1f)});

  // query path that crosses it vertically
  hmap::Path path(
      {hmap::Point(0.5f, 0.1f, 0.8f), hmap::Point(0.5f, 0.9f, 0.2f)});

  hmap::Path truncated = hmap::truncate_at_first_intersection(path, {existing});

  ASSERT_EQ(truncated.size(), 2u);
  EXPECT_NEAR(truncated[1].x, 0.5f, 1e-4f);
  EXPECT_NEAR(truncated[1].y, 0.5f, 1e-4f);
}

TEST(FlowStream, HydrologyPathsToGraph)
{
  // two intersecting streams
  hmap::Path stream1({hmap::Point(0.1f, 0.5f, 1.0f),
                      hmap::Point(0.5f, 0.5f, 0.5f),
                      hmap::Point(0.9f, 0.5f, 0.1f)});

  // stream 2 joins stream 1 at (0.5, 0.5)
  hmap::Path stream2(
      {hmap::Point(0.5f, 0.9f, 0.8f), hmap::Point(0.5f, 0.5f, 0.5f)});

  hmap::Graph g = hmap::hydrology_paths_to_graph({stream1, stream2}, 1e-3f);

  // 4 unique vertices: (0.1, 0.5), (0.5, 0.5), (0.9, 0.5), (0.5, 0.9)
  EXPECT_EQ(g.size(), 4u);
  EXPECT_EQ(g.num_edges(), 3u);

  // compute Strahler order on the resulting graph
  std::vector<int> orders = g.compute_strahler_order();
  ASSERT_EQ(orders.size(), 4u);

  // outlet has minimum elevation (0.1) and should receive the combined flow
  int max_order = *std::max_element(orders.begin(), orders.end());
  EXPECT_EQ(max_order, 2);
}
