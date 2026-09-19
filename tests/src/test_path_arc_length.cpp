#include <cmath>

#include "highmap.hpp"

#include <gtest/gtest.h>

// Open, non-uniformly spaced path: segments of length 1 then 3.
TEST(PathArcLength, OpenNonUniformCumulativeDistance)
{
  hmap::Path path(std::vector<hmap::Point>{{0.f, 0.f, 0.f},
                                           {1.f, 0.f, 0.f},
                                           {4.f, 0.f, 0.f}});

  std::vector<float> cdist = path.get_cumulative_distance();

  ASSERT_EQ(cdist.size(), 3u);
  EXPECT_NEAR(cdist[0], 0.f, 1e-5f);
  EXPECT_NEAR(cdist[1], 1.f, 1e-5f); // first segment length
  EXPECT_NEAR(cdist[2], 4.f, 1e-5f); // total length

  std::vector<float> arc = path.get_arc_length();
  ASSERT_EQ(arc.size(), 3u);
  EXPECT_NEAR(arc[0], 0.f, 1e-5f);
  EXPECT_NEAR(arc[1], 0.25f, 1e-5f);
  EXPECT_NEAR(arc[2], 1.f, 1e-5f);
}

// Two-point open path: must not collapse to a zero-length (NaN) arc.
TEST(PathArcLength, OpenTwoPointArcLengthIsFinite)
{
  hmap::Path path(
      std::vector<hmap::Point>{{0.2f, 0.5f, 0.f}, {0.8f, 0.5f, 0.f}});

  std::vector<float> cdist = path.get_cumulative_distance();
  ASSERT_EQ(cdist.size(), 2u);
  EXPECT_NEAR(cdist[0], 0.f, 1e-5f);
  EXPECT_NEAR(cdist[1], 0.6f, 1e-5f);

  std::vector<float> arc = path.get_arc_length();
  ASSERT_EQ(arc.size(), 2u);
  EXPECT_TRUE(std::isfinite(arc[0]));
  EXPECT_TRUE(std::isfinite(arc[1]));
  EXPECT_NEAR(arc[0], 0.f, 1e-5f);
  EXPECT_NEAR(arc[1], 1.f, 1e-5f);
}

// Closed unit square: cumulative distance must follow consecutive edges and
// include the closing segment (size N + 1).
TEST(PathArcLength, ClosedSquareCumulativeDistance)
{
  hmap::Path path(std::vector<hmap::Point>{{0.f, 0.f, 0.f},
                                           {1.f, 0.f, 0.f},
                                           {1.f, 1.f, 0.f},
                                           {0.f, 1.f, 0.f}});
  path.set_closed(true);

  std::vector<float> cdist = path.get_cumulative_distance();

  ASSERT_EQ(cdist.size(), 5u);
  EXPECT_NEAR(cdist[0], 0.f, 1e-5f);
  EXPECT_NEAR(cdist[1], 1.f, 1e-5f);
  EXPECT_NEAR(cdist[2], 2.f, 1e-5f);
  EXPECT_NEAR(cdist[3], 3.f, 1e-5f);
  EXPECT_NEAR(cdist[4], 4.f, 1e-5f); // full perimeter
}

// Zero-length / degenerate path: arc length must return all zeros without NaNs.
TEST(PathArcLength, DegenerateZeroLengthPath)
{
  hmap::Path path(std::vector<hmap::Point>{{1.f, 1.f, 0.f},
                                           {1.f, 1.f, 0.f},
                                           {1.f, 1.f, 0.f}});

  std::vector<float> arc = path.get_arc_length();
  ASSERT_EQ(arc.size(), 3u);
  EXPECT_FALSE(std::isnan(arc[0]));
  EXPECT_FALSE(std::isnan(arc[1]));
  EXPECT_FALSE(std::isnan(arc[2]));
  EXPECT_NEAR(arc[0], 0.f, 1e-5f);
  EXPECT_NEAR(arc[1], 0.f, 1e-5f);
  EXPECT_NEAR(arc[2], 0.f, 1e-5f);
}

// Empty path edge cases
TEST(PathArcLength, EmptyPathSafety)
{
  hmap::Path path;

  EXPECT_TRUE(path.get_arc_length().empty());
  EXPECT_TRUE(path.get_cumulative_distance().empty());
  EXPECT_TRUE(path.get_curvature().empty());
  EXPECT_TRUE(path.get_tangents().empty());
  EXPECT_TRUE(path.get_normals().empty());
}

// Tangents and Curvature on short paths
TEST(PathArcLength, ShortPathsTangentsAndCurvature)
{
  hmap::Path single_pt(std::vector<hmap::Point>{{1.f, 2.f, 3.f}});
  auto       tg1 = single_pt.get_tangents();
  ASSERT_EQ(tg1.size(), 1u);
  EXPECT_NEAR(tg1[0].x, 1.f, 1e-5f);

  hmap::Path two_pts(
      std::vector<hmap::Point>{{0.f, 0.f, 0.f}, {10.f, 0.f, 0.f}});
  auto tg2 = two_pts.get_tangents();
  ASSERT_EQ(tg2.size(), 2u);
  EXPECT_NEAR(tg2[0].x, 1.f, 1e-5f);
  EXPECT_NEAR(tg2[0].y, 0.f, 1e-5f);
  EXPECT_NEAR(tg2[1].x, 1.f, 1e-5f);
  EXPECT_NEAR(tg2[1].y, 0.f, 1e-5f);

  auto cv2 = two_pts.get_curvature(true);
  ASSERT_EQ(cv2.size(), 2u);
  EXPECT_FALSE(std::isnan(cv2[0]));
  EXPECT_FALSE(std::isnan(cv2[1]));
}
