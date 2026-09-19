#include <vector>

#include <glm/glm.hpp>

#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/kd_tree.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------

static constexpr float eps = 1e-4f;

static bool float_eq(float a, float b, float tol = eps)
{
  return std::abs(a - b) < tol;
}

// ------------------------------------------------------------
// Empty & Single Point
// ------------------------------------------------------------

TEST(KDTreeTest, EmptyTree)
{
  KDTree tree;
  EXPECT_TRUE(tree.empty());
  EXPECT_EQ(tree.size(), 0u);

  EXPECT_EQ(tree.nearest(0.f, 0.f), 0u);

  auto [idx, d_sq] = tree.nearest_with_distance_squared(1.f, 2.f);
  EXPECT_EQ(idx, 0u);
  EXPECT_TRUE(float_eq(d_sq, 0.f));

  std::vector<size_t> indices;
  std::vector<float>  distances;
  tree.neighbor_search(0.f, 0.f, 5, indices, distances);
  EXPECT_TRUE(indices.empty());
  EXPECT_TRUE(distances.empty());

  auto matches = tree.radius_search(0.f, 0.f, 1.f);
  EXPECT_TRUE(matches.empty());

  glm::vec2 drange = tree.compute_neighbor_distance_range(4);
  EXPECT_TRUE(float_eq(drange.x, 0.f));
  EXPECT_TRUE(float_eq(drange.y, 0.f));
}

TEST(KDTreeTest, SinglePoint)
{
  std::vector<float> x = {5.f};
  std::vector<float> y = {10.f};

  KDTree tree(x, y);
  EXPECT_FALSE(tree.empty());
  EXPECT_EQ(tree.size(), 1u);

  EXPECT_EQ(tree.nearest(5.1f, 10.1f), 0u);

  auto [idx, d_sq] = tree.nearest_with_distance_squared(5.f, 13.f);
  EXPECT_EQ(idx, 0u);
  EXPECT_TRUE(float_eq(d_sq, 9.f)); // (13-10)^2 = 9

  auto matches = tree.radius_search(5.f, 10.f, 1.f);
  ASSERT_EQ(matches.size(), 1u);
  EXPECT_EQ(matches[0].first, 0u);
  EXPECT_TRUE(float_eq(matches[0].second, 0.f));

  // negative radius should return empty
  auto no_matches = tree.radius_search(5.f, 10.f, -1.f);
  EXPECT_TRUE(no_matches.empty());
}

// ------------------------------------------------------------
// Multi-Point Queries
// ------------------------------------------------------------

TEST(KDTreeTest, XYVectorQueries)
{
  std::vector<float> x = {0.f, 1.f, 2.f, 0.f};
  std::vector<float> y = {0.f, 0.f, 0.f, 2.f};

  KDTree tree(x, y);
  EXPECT_EQ(tree.size(), 4u);

  // query nearest to (0.1, 0.1) -> should be index 0 (0,0)
  EXPECT_EQ(tree.nearest(0.1f, 0.1f), 0u);

  // query nearest to (1.9, 0.1) -> should be index 2 (2,0)
  EXPECT_EQ(tree.nearest(1.9f, 0.1f), 2u);

  // query nearest to (0.1, 1.9) -> should be index 3 (0,2)
  EXPECT_EQ(tree.nearest(0.1f, 1.9f), 3u);

  // k-NN search with k=2 from (0.9, 0.1) -> index 1 (1,0), then index 0 or 2
  std::vector<size_t> indices;
  std::vector<float>  distances;
  tree.neighbor_search(0.9f, 0.0f, 2, indices, distances);
  ASSERT_EQ(indices.size(), 2u);
  EXPECT_EQ(indices[0], 1u);
  EXPECT_NEAR(distances[0], 0.01f, eps);

  // radius search around (0,0) with radius 1.5 -> points at (0,0) and (1,0)
  auto matches = tree.radius_search(0.f, 0.f, 1.5f);
  EXPECT_EQ(matches.size(), 2u);
}

TEST(KDTreeTest, PointVectorQueries)
{
  std::vector<Point> pts = {
      Point(0.f, 0.f, 10.f),
      Point(1.f, 0.f, 20.f),
      Point(0.f, 1.f, 30.f),
      Point(1.f, 1.f, 40.f),
  };

  KDTree tree(pts);
  EXPECT_EQ(tree.size(), 4u);

  Point q(0.9f, 0.95f, 0.f);
  EXPECT_EQ(tree.nearest(q), 3u);

  glm::vec2 q_vec(0.05f, 0.95f);
  EXPECT_EQ(tree.nearest(q_vec), 2u);

  auto [idx, dist_sq] = tree.nearest_with_distance_squared(q);
  EXPECT_EQ(idx, 3u);
  EXPECT_NEAR(dist_sq, 0.01f + 0.0025f, eps);

  std::vector<size_t> indices;
  std::vector<float>  distances;
  tree.neighbor_search(q, 4, indices, distances);
  EXPECT_EQ(indices.size(), 4u);
  EXPECT_EQ(indices[0], 3u);
}

// ------------------------------------------------------------
// Cloud and Path Construction
// ------------------------------------------------------------

TEST(KDTreeTest, CloudConstruction)
{
  Cloud cloud;
  cloud.emplace_back(0.f, 0.f);
  cloud.emplace_back(2.f, 2.f);
  cloud.emplace_back(4.f, 4.f);

  KDTree tree(cloud);
  EXPECT_EQ(tree.size(), 3u);
  EXPECT_EQ(tree.nearest(1.9f, 2.1f), 1u);
}

TEST(KDTreeTest, PathConstruction)
{
  std::vector<Point> pts = {Point(0, 0, 0), Point(5, 5, 0), Point(10, 10, 0)};
  Path               path(pts);

  KDTree tree(path);
  EXPECT_EQ(tree.size(), 3u);
  EXPECT_EQ(tree.nearest(4.8f, 5.2f), 1u);
}

// ------------------------------------------------------------
// Distance Range
// ------------------------------------------------------------

TEST(KDTreeTest, DistanceRangeEuclidean)
{
  // 4 corners of unit square: (0,0), (1,0), (0,1), (1,1)
  std::vector<Point> pts = {
      Point(0.f, 0.f, 0.f),
      Point(1.f, 0.f, 0.f),
      Point(0.f, 1.f, 0.f),
      Point(1.f, 1.f, 0.f),
  };

  KDTree tree(pts);

  // For each corner with k=2:
  // neighbor 0 is the point itself (dist = 0)
  // neighbor 1 is distance 1.0
  // So min distance > 0 is 1.0, max distance is 1.0
  glm::vec2 range = tree.compute_neighbor_distance_range(2);
  EXPECT_NEAR(range.x, 1.0f, eps);
  EXPECT_NEAR(range.y, 1.0f, eps);

  // With k=4, the farthest neighbor is at diagonal dist sqrt(2) ~ 1.4142
  glm::vec2 range_diag = tree.compute_neighbor_distance_range(4);
  EXPECT_NEAR(range_diag.x, 1.0f, eps);
  EXPECT_NEAR(range_diag.y, std::sqrt(2.0f), eps);
}

// ------------------------------------------------------------
// Move Semantics
// ------------------------------------------------------------

TEST(KDTreeTest, MoveSemantics)
{
  std::vector<Point> pts = {Point(1.f, 2.f, 0.f), Point(3.f, 4.f, 0.f)};

  KDTree tree1(pts);
  EXPECT_EQ(tree1.size(), 2u);

  KDTree tree2(std::move(tree1));
  EXPECT_EQ(tree2.size(), 2u);
  EXPECT_EQ(tree2.nearest(1.1f, 1.9f), 0u);

  KDTree tree3;
  tree3 = std::move(tree2);
  EXPECT_EQ(tree3.size(), 2u);
  EXPECT_EQ(tree3.nearest(3.1f, 3.9f), 1u);
}
