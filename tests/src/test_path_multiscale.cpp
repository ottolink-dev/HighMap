#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(PathMultiscale, ConnectsStartToEnd)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);

  glm::ivec2 start = {10, 10};
  glm::ivec2 end = {115, 115};

  std::vector<int> ip, jp;
  hmap::find_path_multiscale(z,
                             start,
                             end,
                             ip,
                             jp,
                             3,
                             4,
                             1.f,
                             0.1f,
                             2.f,
                             1.f,
                             nullptr,
                             true);

  ASSERT_FALSE(ip.empty());
  ASSERT_EQ(ip.size(), jp.size());

  // check start and end points
  EXPECT_EQ(ip.front(), start.x);
  EXPECT_EQ(jp.front(), start.y);
  EXPECT_EQ(ip.back(), end.x);
  EXPECT_EQ(jp.back(), end.y);

  // check 8-connectivity
  for (size_t k = 1; k < ip.size(); ++k)
  {
    int di = std::abs(ip[k] - ip[k - 1]);
    int dj = std::abs(jp[k] - jp[k - 1]);
    EXPECT_LE(di, 1);
    EXPECT_LE(dj, 1);
    EXPECT_GT(di + dj, 0);
  }
}

TEST(PathMultiscale, VectorOverloadReturnsSameEndpoints)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z(shape, 0.5f);

  glm::ivec2 start = {5, 5};
  glm::ivec2 end = {58, 58};

  std::vector<glm::ivec2> path = hmap::find_path_multiscale(z,
                                                            start,
                                                            end,
                                                            3,
                                                            3);

  ASSERT_FALSE(path.empty());
  EXPECT_EQ(path.front(), start);
  EXPECT_EQ(path.back(), end);
}

TEST(PathMultiscale, IdenticalStartAndEnd)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z(shape, 1.f);

  glm::ivec2       start = {20, 20};
  std::vector<int> ip, jp;
  hmap::find_path_multiscale(z, start, start, ip, jp);

  ASSERT_EQ(ip.size(), 1u);
  ASSERT_EQ(jp.size(), 1u);
  EXPECT_EQ(ip[0], start.x);
  EXPECT_EQ(jp[0], start.y);
}

TEST(PathMultiscale, AvoidsObstaclesWithNogoMask)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z(shape, 0.f);

  // create a vertical wall in the middle with a single opening at the top
  hmap::Array mask(shape, 0.f);
  for (int j = 0; j < 50; ++j)
  {
    mask(32, j) = 1.f;
  }

  glm::ivec2 start = {10, 25};
  glm::ivec2 end = {50, 25};

  std::vector<glm::ivec2> path = hmap::find_path_multiscale(z,
                                                            start,
                                                            end,
                                                            2,
                                                            8,
                                                            1.f,
                                                            0.1f,
                                                            2.f,
                                                            1.f,
                                                            &mask,
                                                            true);

  ASSERT_FALSE(path.empty());
  EXPECT_EQ(path.front(), start);
  EXPECT_EQ(path.back(), end);

  // path must not cross the wall
  for (const auto &p : path)
  {
    if (p.x == 32)
    {
      EXPECT_GE(p.y, 50);
    }
  }
}

TEST(PathMultiscale, DijkstraVsAStar)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {3.f, 3.f},
                                  123);

  glm::ivec2 start = {15, 15};
  glm::ivec2 end = {110, 110};

  // test both astar and dijkstra modes
  std::vector<glm::ivec2> path_astar = hmap::find_path_multiscale(z,
                                                                  start,
                                                                  end,
                                                                  3,
                                                                  4,
                                                                  1.f,
                                                                  0.1f,
                                                                  2.f,
                                                                  1.f,
                                                                  nullptr,
                                                                  true);
  std::vector<glm::ivec2> path_dijkstra = hmap::find_path_multiscale(z,
                                                                     start,
                                                                     end,
                                                                     3,
                                                                     4,
                                                                     1.f,
                                                                     0.1f,
                                                                     2.f,
                                                                     1.f,
                                                                     nullptr,
                                                                     false);

  ASSERT_FALSE(path_astar.empty());
  ASSERT_FALSE(path_dijkstra.empty());

  EXPECT_EQ(path_astar.front(), start);
  EXPECT_EQ(path_astar.back(), end);
  EXPECT_EQ(path_dijkstra.front(), start);
  EXPECT_EQ(path_dijkstra.back(), end);
}

TEST(PathMultiscale, OptionalPathSmoothing)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  7);

  glm::ivec2 start = {10, 10};
  glm::ivec2 end = {110, 110};

  std::vector<glm::ivec2> path_smooth = hmap::find_path_multiscale(z,
                                                                   start,
                                                                   end,
                                                                   3,
                                                                   4,
                                                                   1.f,
                                                                   0.1f,
                                                                   2.f,
                                                                   1.f,
                                                                   nullptr,
                                                                   true,
                                                                   true);

  ASSERT_FALSE(path_smooth.empty());
  EXPECT_EQ(path_smooth.front(), start);
  EXPECT_EQ(path_smooth.back(), end);
}

TEST(PathMultiscale, PathOverloadAndToArray)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  99);

  glm::ivec2 start = {8, 8};
  glm::ivec2 end = {55, 55};
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  hmap::Path path = hmap::find_path_multiscale(z, start, end, bbox, 3, 4);

  ASSERT_FALSE(path.points.empty());

  // verify point values hold actual elevations from z
  for (const auto &p : path.points)
  {
    int i = int(std::round(p.x * float(shape.x - 1)));
    int j = int(std::round(p.y * float(shape.y - 1)));
    i = std::clamp(i, 0, shape.x - 1);
    j = std::clamp(j, 0, shape.y - 1);
    EXPECT_NEAR(p.v, z(i, j), 1e-4f);
  }

  // verify exact elevation on path endpoints
  EXPECT_NEAR(path.points.front().v, z(start.x, start.y), 1e-5f);
  EXPECT_NEAR(path.points.back().v, z(end.x, end.y), 1e-5f);

  // test built-in to_array
  hmap::Array arr = path.to_array(shape, bbox);
  EXPECT_EQ(arr.shape, shape);

  // verify start and end cells are rasterized in exported array
  EXPECT_NEAR(arr(start.x, start.y), z(start.x, start.y), 0.05f);
  EXPECT_NEAR(arr(end.x, end.y), z(end.x, end.y), 0.05f);
}

TEST(PathMultiscale, FindCutPathMultiscale)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);

  hmap::Path path = hmap::find_cut_path_multiscale(
      z,
      hmap::DomainBoundary::BOUNDARY_LEFT,
      hmap::DomainBoundary::BOUNDARY_RIGHT,
      3,
      6);

  ASSERT_FALSE(path.points.empty());
  EXPECT_NEAR(path.points.front().x, 0.f, 1e-4f);
  EXPECT_NEAR(path.points.back().x, 1.f, 1e-4f);
}

TEST(PathMultiscale, FindPathDijkstraHarmonized)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  10);

  glm::ivec2 start = {5, 5};
  glm::ivec2 end = {58, 58};
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  // vector overload
  std::vector<glm::ivec2> vec_path = hmap::find_path_dijkstra(z, start, end);
  ASSERT_FALSE(vec_path.empty());
  EXPECT_EQ(vec_path.front(), start);
  EXPECT_EQ(vec_path.back(), end);

  // Path overload
  hmap::Path path = hmap::find_path_dijkstra(z, start, end, bbox);
  ASSERT_FALSE(path.points.empty());
  EXPECT_NEAR(path.points.front().v, z(start.x, start.y), 1e-5f);
  EXPECT_NEAR(path.points.back().v, z(end.x, end.y), 1e-5f);
}

TEST(PathMultiscale, FindPathMidpointHarmonized)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  20);

  glm::ivec2 start = {5, 5};
  glm::ivec2 end = {58, 58};
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  // void (i_path, j_path) overload
  std::vector<int> ip, jp;
  hmap::find_path_midpoint(z, start, end, ip, jp);
  ASSERT_FALSE(ip.empty());
  EXPECT_EQ(ip.front(), start.x);
  EXPECT_EQ(jp.front(), start.y);
  EXPECT_EQ(ip.back(), end.x);
  EXPECT_EQ(jp.back(), end.y);

  // Path overload
  hmap::Path path = hmap::find_path_midpoint(z, start, end, bbox);
  ASSERT_FALSE(path.points.empty());
  EXPECT_NEAR(path.points.front().v, z(start.x, start.y), 1e-5f);
  EXPECT_NEAR(path.points.back().v, z(end.x, end.y), 1e-5f);
}
