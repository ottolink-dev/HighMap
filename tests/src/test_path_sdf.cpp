#include <cmath>
#include <vector>

#include "highmap/geometry/path.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(PathSDF, OpenPathDoesNotIncludeClosingSegment)
{
  Path path(std::vector<Point>{{0.25f, 0.25f, 0.f},
                               {0.75f, 0.25f, 0.f},
                               {0.75f, 0.75f, 0.f},
                               {0.25f, 0.75f, 0.f}});
  path.set_closed(false);

  Array sdf = path_sdf_to_array(path, {8, 8});

  // The missing left edge must remain a gap, with unsigned distances.
  EXPECT_NEAR(sdf(2, 4), 0.25f, 1e-6f);
  EXPECT_NEAR(sdf(1, 4), std::sqrt(0.125f * 0.125f + 0.25f * 0.25f), 1e-6f);
  EXPECT_NEAR(sdf(4, 4), 0.25f, 1e-6f);
  EXPECT_NEAR(sdf(4, 2), 0.f, 1e-6f);
  EXPECT_NEAR(sdf(6, 4), 0.f, 1e-6f);
  EXPECT_NEAR(sdf(4, 6), 0.f, 1e-6f);
}

TEST(PathSDF, ClosedPathIncludesClosingSegmentAndInteriorSign)
{
  Path path(std::vector<Point>{{0.25f, 0.25f, 0.f},
                               {0.75f, 0.25f, 0.f},
                               {0.75f, 0.75f, 0.f},
                               {0.25f, 0.75f, 0.f}});
  path.set_closed(true);

  Array sdf = path_sdf_to_array(path, {8, 8});

  EXPECT_NEAR(sdf(2, 4), 0.f, 1e-6f);
  EXPECT_NEAR(sdf(1, 4), 0.125f, 1e-6f);
  EXPECT_NEAR(sdf(4, 4), -0.25f, 1e-6f);
}

TEST(PathSDF, TwoPointOpenPathIncludesSegmentAndEndpoints)
{
  Path path(std::vector<Point>{{0.25f, 0.5f, 0.f}, {0.75f, 0.5f, 0.f}});
  path.set_closed(false);

  Array sdf = path_sdf_to_array(path, {8, 8});

  EXPECT_NEAR(sdf(4, 4), 0.f, 1e-6f);
  EXPECT_NEAR(sdf(4, 2), 0.25f, 1e-6f);
  EXPECT_NEAR(sdf(1, 4), 0.125f, 1e-6f);
  EXPECT_NEAR(sdf(7, 4), 0.125f, 1e-6f);
}

TEST(PathSDF, OpenPathRespectsDomainAndWarp)
{
  Path path(std::vector<Point>{{2.f, 3.f, 0.f},
                               {4.f, 3.f, 0.f},
                               {4.f, 5.f, 0.f},
                               {2.f, 5.f, 0.f}});
  path.set_closed(false);
  Array dx({8, 8}, 0.25f);
  Array dy({8, 8}, -0.25f);

  Array sdf = path_sdf_to_array(path,
                                {8, 8},
                                BBox2D({1.f, 5.f, 2.f, 6.f}),
                                &dx,
                                &dy);

  // Grid point (1.5, 4) is warped to (1.75, 3.75), nearest the first endpoint.
  EXPECT_NEAR(sdf(1, 4), std::sqrt(0.25f * 0.25f + 0.75f * 0.75f), 1e-6f);
}
