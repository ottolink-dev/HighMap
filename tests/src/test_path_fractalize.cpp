#include "highmap.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(PathFractalize, BoundedWithinSingleEdgeBoundingBox)
{
  Point p1(0.f, 0.f);
  Point p2(10.f, 10.f);
  Path  path({p1, p2});

  int           iterations = 5;
  std::uint32_t seed = 42;
  float         sigma = 0.5f;

  Path result_bounded = fractalize(path,
                                   iterations,
                                   seed,
                                   sigma,
                                   0,
                                   1.f,
                                   nullptr,
                                   {0.f, 10.f, 0.f, 10.f},
                                   true);

  for (const auto &pt : result_bounded.points)
  {
    EXPECT_GE(pt.x, 0.f);
    EXPECT_LE(pt.x, 10.f);
    EXPECT_GE(pt.y, 0.f);
    EXPECT_LE(pt.y, 10.f);
  }
}

TEST(PathFractalize, BoundedMultiEdgePath)
{
  // Edge 0: (0, 0) -> (5, 10), bbox = [0, 5] x [0, 10]
  // Edge 1: (5, 10) -> (10, 0), bbox = [5, 10] x [0, 10]
  Path path({Point(0.f, 0.f), Point(5.f, 10.f), Point(10.f, 0.f)});

  int           iterations = 4;
  std::uint32_t seed = 99;
  float         sigma = 0.5f;

  Path result = fractalize(path,
                           iterations,
                           seed,
                           sigma,
                           0,
                           1.f,
                           nullptr,
                           {0.f, 10.f, 0.f, 10.f},
                           true);

  // At iteration `it`, each initial edge is subdivided into 2^it segments,
  // having (2^it + 1) points. Total points for 2 open edges = 2 * 2^it + 1 = 2
  // * 16 + 1 = 33 points.
  size_t points_per_edge = size_t(1) << iterations;

  for (size_t i = 0; i <= points_per_edge; ++i)
  {
    const auto &pt = result.points[i];
    EXPECT_GE(pt.x, 0.f);
    EXPECT_LE(pt.x, 5.f);
    EXPECT_GE(pt.y, 0.f);
    EXPECT_LE(pt.y, 10.f);
  }

  for (size_t i = points_per_edge; i < result.points.size(); ++i)
  {
    const auto &pt = result.points[i];
    EXPECT_GE(pt.x, 5.f);
    EXPECT_LE(pt.x, 10.f);
    EXPECT_GE(pt.y, 0.f);
    EXPECT_LE(pt.y, 10.f);
  }
}

TEST(PathFractalize, DefaultUnboundedBehavior)
{
  Point p1(0.f, 0.f);
  Point p2(10.f, 0.f);
  Path  path({p1, p2});

  int           iterations = 3;
  std::uint32_t seed = 123;
  float         sigma = 0.5f;

  // Unbounded: perpendicular displacement will generate points with non-zero y
  Path result_unbounded = fractalize(path, iterations, seed, sigma);

  bool has_non_zero_y = false;
  for (const auto &pt : result_unbounded.points)
  {
    if (std::abs(pt.y) > 1e-4f)
    {
      has_non_zero_y = true;
      break;
    }
  }
  EXPECT_TRUE(has_non_zero_y);

  // Bounded: clamped to min(p1.y, p2.y) = 0 and max(p1.y, p2.y) = 0, so all y
  // must remain 0
  Path result_bounded = fractalize(path,
                                   iterations,
                                   seed,
                                   sigma,
                                   0,
                                   1.f,
                                   nullptr,
                                   {0.f, 10.f, 0.f, 0.f},
                                   true);

  for (const auto &pt : result_bounded.points)
  {
    EXPECT_NEAR(pt.y, 0.f, 1e-5f);
    EXPECT_GE(pt.x, 0.f);
    EXPECT_LE(pt.x, 10.f);
  }
}

TEST(PathFractalizeUniform, UniformDisplacementAlongDifferentLengthEdges)
{
  // Edge 1: length = 1.0 (from (0,0) to (1,0))
  // Edge 2: length = 10.0 (from (1,0) to (11,0))
  Path path({Point(0.f, 0.f), Point(1.f, 0.f), Point(11.f, 0.f)});

  float sigma = 0.2f;
  float spacing = 0.5f;
  Path  res = fractalize_uniform(path, 3, 42, sigma, spacing);

  EXPECT_GT(res.size(), path.size());

  for (const auto &pt : res.points)
  {
    EXPECT_GE(pt.x, -1.f);
    EXPECT_LE(pt.x, 12.f);
  }
}
