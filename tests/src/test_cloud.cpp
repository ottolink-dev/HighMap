#include <glm/glm.hpp>

#include "highmap/geometry/cloud.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------

static constexpr float eps = 1e-5f;

bool floatEq(float a, float b, float tol = eps)
{
  return std::abs(a - b) < tol;
}

void expectVecNear(const std::vector<float> &a,
                   const std::vector<float> &b,
                   float                     tol = eps)
{
  ASSERT_EQ(a.size(), b.size());
  for (size_t i = 0; i < a.size(); ++i)
    EXPECT_NEAR(a[i], b[i], tol);
}

// ------------------------------------------------------------
// Constructors
// ------------------------------------------------------------

TEST(CloudTest, DefaultConstructorEmpty)
{
  Cloud cloud;
  EXPECT_EQ(cloud.size(), 0);
}

TEST(CloudTest, ConstructFromXY)
{
  std::vector<float> x = {0.f, 1.f};
  std::vector<float> y = {2.f, 3.f};

  Cloud cloud(x, y);

  ASSERT_EQ(cloud.size(), 2);
  EXPECT_TRUE(floatEq(cloud[0].x, 0.f));
  EXPECT_TRUE(floatEq(cloud[0].y, 2.f));
}

TEST(CloudTest, ConstructFromXYV)
{
  std::vector<float> x = {0.f, 1.f};
  std::vector<float> y = {2.f, 3.f};
  std::vector<float> v = {10.f, 20.f};

  Cloud cloud(x, y, v);

  ASSERT_EQ(cloud.size(), 2);
  EXPECT_TRUE(floatEq(cloud[1].v, 20.f));
}

TEST(CloudTest, RandomConstructorDeterministic)
{
  Cloud c1(10, 42);
  Cloud c2(10, 42);

  ASSERT_EQ(c1.size(), c2.size());

  for (size_t i = 0; i < c1.size(); ++i)
  {
    EXPECT_NEAR(c1[i].x, c2[i].x, eps);
    EXPECT_NEAR(c1[i].y, c2[i].y, eps);
  }
}

// ------------------------------------------------------------
// Basic Operations
// ------------------------------------------------------------

TEST(CloudTest, PushBackIncreasesSize)
{
  Cloud cloud;

  cloud.push_back({1.f, 2.f, 3.f});

  EXPECT_EQ(cloud.size(), 1);
}

TEST(CloudTest, EraseReducesSize)
{
  Cloud cloud(std::vector<Point>{{0, 0, 0}, {1, 1, 1}});

  cloud.erase(cloud.begin());

  ASSERT_EQ(cloud.size(), 1);
  EXPECT_TRUE(floatEq(cloud[0].x, 1.f));
}

TEST(CloudTest, ClearEmptiesCloud)
{
  Cloud cloud(std::vector<Point>{{0, 0, 0}, {1, 1, 1}});

  cloud.clear();

  EXPECT_EQ(cloud.size(), 0);
}

// ------------------------------------------------------------
// Accessors
// ------------------------------------------------------------

TEST(CloudTest, GetBBoxCorrect)
{
  Cloud cloud(std::vector<Point>{{-1.f, 2.f, 0.f}, {3.f, 4.f, 0.f}});

  glm::vec4 bbox = cloud.get_bbox();

  EXPECT_TRUE(floatEq(bbox.x, -1.f));
  EXPECT_TRUE(floatEq(bbox.y, 3.f));
  EXPECT_TRUE(floatEq(bbox.z, 2.f));
  EXPECT_TRUE(floatEq(bbox.w, 4.f));
}

TEST(CloudTest, GetCenterCorrect)
{
  Cloud cloud(std::vector<Point>{{0.f, 0.f, 0.f}, {2.f, 2.f, 0.f}});

  Point c = cloud.get_center();

  EXPECT_TRUE(floatEq(c.x, 1.f));
  EXPECT_TRUE(floatEq(c.y, 1.f));
}

TEST(CloudTest, GetValuesCorrect)
{
  Cloud cloud(std::vector<Point>{{0, 0, 1}, {1, 1, 2}});

  auto values = cloud.get_values();

  expectVecNear(values, {1.f, 2.f});
}

TEST(CloudTest, GetValuesMinMax)
{
  Cloud cloud(std::vector<Point>{{0, 0, 5}, {1, 1, 2}, {2, 2, 9}});

  EXPECT_TRUE(floatEq(cloud.get_values_min(), 2.f));
  EXPECT_TRUE(floatEq(cloud.get_values_max(), 9.f));
}

TEST(CloudTest, GetXYConcatenation)
{
  Cloud cloud(std::vector<Point>{{1, 2, 0}, {3, 4, 0}});

  auto xy = cloud.get_xy();

  expectVecNear(xy, {1, 2, 3, 4});
}

// ------------------------------------------------------------
// Nearest Point
// ------------------------------------------------------------

TEST(CloudTest, NearestPointCorrect)
{
  Cloud cloud(std::vector<Point>{{0, 0, 0}, {10, 10, 0}});

  size_t idx = cloud.nearest_point({1.f, 1.f});

  EXPECT_EQ(idx, 0);
}

// ------------------------------------------------------------
// Value Operations
// ------------------------------------------------------------

TEST(CloudTest, SetValuesConstant)
{
  Cloud cloud(std::vector<Point>{{0, 0, 1}, {1, 1, 2}});

  cloud.set_values(5.f);

  for (auto &p : cloud)
    EXPECT_TRUE(floatEq(p.v, 5.f));
}

TEST(CloudTest, SetValuesVector)
{
  Cloud cloud(std::vector<Point>{{0, 0, 0}, {1, 1, 0}});

  cloud.set_values({3.f, 4.f});

  EXPECT_TRUE(floatEq(cloud[0].v, 3.f));
  EXPECT_TRUE(floatEq(cloud[1].v, 4.f));
}

// ------------------------------------------------------------
// Remap
// ------------------------------------------------------------

TEST(CloudTest, RemapValuesRange)
{
  Cloud cloud(std::vector<Point>{{0, 0, 0}, {1, 1, 10}});

  cloud.remap_values(0.f, 1.f);

  EXPECT_TRUE(floatEq(cloud.get_values_min(), 0.f));
  EXPECT_TRUE(floatEq(cloud.get_values_max(), 1.f));
}

// ------------------------------------------------------------
// Randomization / Shuffle
// ------------------------------------------------------------

TEST(CloudTest, RandomizeDeterministic)
{
  Cloud c1(std::vector<Point>{{0, 0, 0}, {1, 1, 1}});
  Cloud c2 = c1;

  c1.randomize(42);
  c2.randomize(42);

  for (size_t i = 0; i < c1.size(); ++i)
  {
    EXPECT_NEAR(c1[i].x, c2[i].x, eps);
    EXPECT_NEAR(c1[i].y, c2[i].y, eps);
  }
}

TEST(CloudTest, ShuffleDeterministic)
{
  Cloud c1(std::vector<Point>{{0, 0, 0}, {1, 1, 1}});
  Cloud c2 = c1;

  c1.shuffle(1.f, 1.f, 42, 1.f);
  c2.shuffle(1.f, 1.f, 42, 1.f);

  for (size_t i = 0; i < c1.size(); ++i)
  {
    EXPECT_NEAR(c1[i].x, c2[i].x, eps);
    EXPECT_NEAR(c1[i].v, c2[i].v, eps);
  }
}

// ------------------------------------------------------------
// Convex Hull
// ------------------------------------------------------------

TEST(CloudTest, ConvexHullSquare)
{
  Cloud cloud(std::vector<Point>{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}});

  auto hull = cloud.get_convex_hull();

  EXPECT_EQ(hull.size(), 4);
}

TEST(CloudTest, ConvexHullIgnoresInnerPoint)
{
  Cloud cloud(std::vector<Point>{{0, 0, 0},
                                 {1, 0, 0},
                                 {1, 1, 0},
                                 {0, 1, 0},
                                 {0.5f, 0.5f, 0}});

  auto hull = cloud.get_convex_hull();

  EXPECT_EQ(hull.size(), 4);
}

// ------------------------------------------------------------
// Merge (Free functions)
// ------------------------------------------------------------

TEST(CloudTest, MergeCloudsConcatenates)
{
  Cloud a(std::vector<Point>{{0, 0, 0}});
  Cloud b(std::vector<Point>{{1, 1, 1}});

  Cloud merged = merge_cloud(a, b);

  EXPECT_EQ(merged.size(), 2);
}

// ------------------------------------------------------------
// Scale (Free functions)
// ------------------------------------------------------------

TEST(CloudTest, ScaleUniformDefaultCenter)
{
  Cloud cloud(std::vector<Point>{{0.f, 0.f, 10.f}, {1.f, 1.f, 20.f}});

  // Default center of unit bbox is (0.5, 0.5). Scale 0.8 should shrink towards
  // (0.5, 0.5) x: 0.5 + 0.8 * (0 - 0.5) = 0.1 y: 0.5 + 0.8 * (1 - 0.5) = 0.9
  Cloud scaled = scale(cloud, 0.8f);

  ASSERT_EQ(scaled.size(), 2);
  EXPECT_NEAR(scaled[0].x, 0.1f, eps);
  EXPECT_NEAR(scaled[0].y, 0.1f, eps);
  EXPECT_NEAR(scaled[0].v, 10.f, eps);

  EXPECT_NEAR(scaled[1].x, 0.9f, eps);
  EXPECT_NEAR(scaled[1].y, 0.9f, eps);
  EXPECT_NEAR(scaled[1].v, 20.f, eps);
}

TEST(CloudTest, ScaleNonUniformCustomCenter)
{
  Cloud cloud(std::vector<Point>{{1.f, 2.f, 5.f}});

  // Scale x by 2, y by 0.5 relative to (1.f, 1.f)
  // x: 1 + 2 * (1 - 1) = 1
  // y: 1 + 0.5 * (2 - 1) = 1.5
  Cloud scaled = scale(cloud, glm::vec2(2.f, 0.5f), glm::vec2(1.f, 1.f));

  ASSERT_EQ(scaled.size(), 1);
  EXPECT_NEAR(scaled[0].x, 1.f, eps);
  EXPECT_NEAR(scaled[0].y, 1.5f, eps);
  EXPECT_NEAR(scaled[0].v, 5.f, eps);
}

TEST(CloudTest, ScaleCustomCenter)
{
  Cloud     cloud(std::vector<Point>{{10.f, 20.f, 1.f}, {20.f, 40.f, 2.f}});
  glm::vec2 center = {15.f, 30.f};

  // Scale 0.5 relative to custom center (15, 30)
  // p0: x = 15 + 0.5 * (10 - 15) = 12.5, y = 30 + 0.5 * (20 - 30) = 25
  // p1: x = 15 + 0.5 * (20 - 15) = 17.5, y = 30 + 0.5 * (40 - 30) = 35
  Cloud scaled = scale(cloud, 0.5f, center);

  ASSERT_EQ(scaled.size(), 2);
  EXPECT_NEAR(scaled[0].x, 12.5f, eps);
  EXPECT_NEAR(scaled[0].y, 25.f, eps);
  EXPECT_NEAR(scaled[1].x, 17.5f, eps);
  EXPECT_NEAR(scaled[1].y, 35.f, eps);
}

TEST(CloudTest, ScalePreservesPointValues)
{
  Cloud cloud(std::vector<Point>{
      {0.1f, 0.2f, 100.5f},
      {0.5f, 0.5f, -42.0f},
      {0.9f, 0.8f, 3.1415f},
  });

  Cloud scaled = scale(cloud, glm::vec2(0.5f, 1.5f));

  ASSERT_EQ(scaled.size(), 3u);
  EXPECT_FLOAT_EQ(scaled[0].v, 100.5f);
  EXPECT_FLOAT_EQ(scaled[1].v, -42.0f);
  EXPECT_FLOAT_EQ(scaled[2].v, 3.1415f);
}

// ------------------------------------------------------------
// Container Interface & Move Semantics
// ------------------------------------------------------------

TEST(CloudTest, MoveConstructor)
{
  std::vector<Point> pts = {{1.f, 2.f, 3.f}, {4.f, 5.f, 6.f}};
  Cloud              cloud(std::move(pts));

  EXPECT_EQ(cloud.size(), 2);
  EXPECT_TRUE(floatEq(cloud[0].x, 1.f));
  EXPECT_TRUE(floatEq(cloud[1].v, 6.f));
}

TEST(CloudTest, ContainerAccessAndIterators)
{
  Cloud cloud(std::vector<Point>{{1.f, 2.f, 3.f}, {4.f, 5.f, 6.f}});

  // Indexing
  EXPECT_TRUE(floatEq(cloud[0].x, 1.f));
  EXPECT_TRUE(floatEq(cloud.at(1).y, 5.f));
  EXPECT_TRUE(floatEq(cloud.front().v, 3.f));
  EXPECT_TRUE(floatEq(cloud.back().x, 4.f));

  // Mutation via reference
  cloud[0].x = 10.f;
  EXPECT_TRUE(floatEq(cloud[0].x, 10.f));

  // Range-based for loop
  float sum_v = 0.f;
  for (const auto &p : cloud)
    sum_v += p.v;
  EXPECT_TRUE(floatEq(sum_v, 9.f));

  // Capacity, emplace_back & push_back
  cloud.reserve(10);
  EXPECT_GE(cloud.capacity(), 10u);
  cloud.emplace_back(7.f, 8.f, 9.f);
  EXPECT_EQ(cloud.size(), 3u);
  EXPECT_TRUE(floatEq(cloud.back().v, 9.f));

  Point p_lvalue(10.f, 11.f, 12.f);
  cloud.push_back(p_lvalue);
  cloud.push_back(Point(13.f, 14.f, 15.f));
  EXPECT_EQ(cloud.size(), 5u);
  EXPECT_TRUE(floatEq(cloud.back().x, 13.f));

  auto it = cloud.erase(cloud.begin() + 1);
  EXPECT_EQ(cloud.size(), 4u);
  EXPECT_TRUE(floatEq(it->x, 7.f));
}

// ------------------------------------------------------------
// Rejection Filter & Edge Cases
// ------------------------------------------------------------

TEST(CloudTest, RejectionFilterDensityActuallyErases)
{
  // Cloud of 10 points
  std::vector<Point> pts;
  for (int i = 0; i < 10; ++i)
    pts.emplace_back(0.5f, 0.5f, 0.f);
  Cloud cloud(pts);

  // Mask where density is 0 everywhere -> all points should be removed
  Array zero_mask({10, 10}, 0.f);
  rejection_filter_density(cloud, zero_mask, 42);

  EXPECT_EQ(cloud.size(), 0u);
  EXPECT_TRUE(cloud.empty());
}

TEST(CloudTest, RemapValuesEmptyCloudSafe)
{
  Cloud cloud;
  EXPECT_NO_THROW(cloud.remap_values(0.f, 1.f));
}

TEST(CloudTest, ConvexHullDegenerateInputsSafe)
{
  Cloud empty_cloud;
  EXPECT_TRUE(empty_cloud.get_convex_hull().empty());

  Cloud two_points(std::vector<Point>{{0.f, 0.f, 0.f}, {1.f, 1.f, 0.f}});
  EXPECT_TRUE(two_points.get_convex_hull().empty());
}

TEST(CloudTest, ToArrayDegenerateBBoxSafe)
{
  Cloud cloud(std::vector<Point>{{0.5f, 0.5f, 1.f}});
  Array array({10, 10}, 0.f);

  // Degenerate bbox (width = 0, height = 0)
  EXPECT_NO_THROW(cloud.to_array(array, {0.f, 0.f, 0.f, 0.f}));
}
