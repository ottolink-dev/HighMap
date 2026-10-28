#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/sdf.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using namespace hmap;

TEST(BVH2DPathSDF, OpenPathMatchesBruteForceAndBVH)
{
  Path path(std::vector<Point>{{0.2f, 0.3f, 0.f},
                               {0.5f, 0.8f, 0.f},
                               {0.8f, 0.4f, 0.f},
                               {0.9f, 0.7f, 0.f}});
  path.set_closed(false);

  glm::ivec2 shape(32, 32);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array bf = path_sdf_to_array_brute_force(path, shape, bbox);
  Array bvh = path_sdf_to_array_bvh(path, shape, bbox);
  Array def = path_sdf_to_array(path, shape, bbox);
  Array no_bvh = path_sdf_to_array(path, shape, bbox, nullptr, nullptr, false);

  ASSERT_EQ(bf.shape.x, bvh.shape.x);
  ASSERT_EQ(bf.shape.y, bvh.shape.y);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-5f)
          << "Mismatch at pixel (" << i << ", " << j << ")";
      EXPECT_NEAR(def(i, j), bvh(i, j), 1e-5f);
      EXPECT_NEAR(no_bvh(i, j), bf(i, j), 1e-5f);
    }
  }
}

TEST(BVH2DPathSDF, ClosedPathMatchesBruteForceAndBVH)
{
  Path path(std::vector<Point>{{0.25f, 0.25f, 0.f},
                               {0.75f, 0.25f, 0.f},
                               {0.75f, 0.75f, 0.f},
                               {0.25f, 0.75f, 0.f}});
  path.set_closed(true);

  glm::ivec2 shape(24, 24);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array bf = path_sdf_to_array_brute_force(path, shape, bbox);
  Array bvh = path_sdf_to_array_bvh(path, shape, bbox);
  Array def = path_sdf_to_array(path, shape, bbox);

  ASSERT_EQ(bf.shape.x, bvh.shape.x);
  ASSERT_EQ(bf.shape.y, bvh.shape.y);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-5f)
          << "Closed path mismatch at pixel (" << i << ", " << j << ")";
      EXPECT_NEAR(def(i, j), bvh(i, j), 1e-5f);
    }
  }
}

TEST(BVH2DPathSDF, PolylineFunctionsEquivalenceAndWarp)
{
  Path path(std::vector<Point>{{2.f, 3.f, 0.f},
                               {3.5f, 4.2f, 0.f},
                               {4.f, 3.f, 0.f},
                               {4.f, 5.f, 0.f},
                               {2.f, 5.f, 0.f}});
  path.set_closed(false);

  glm::ivec2 shape(24, 24);
  glm::vec4  bbox(1.f, 5.f, 2.f, 6.f);

  Array dx(shape, 0.15f);
  Array dy(shape, -0.1f);

  Array bf = sdf_2d_polyline_brute_force(path, shape, bbox, &dx, &dy);
  Array bvh = sdf_2d_polyline_bvh(path, shape, bbox, &dx, &dy);
  Array def = sdf_2d_polyline(path, shape, bbox, &dx, &dy);
  Array no_bvh = sdf_2d_polyline(path, shape, bbox, &dx, &dy, false);

  ASSERT_EQ(bf.shape.x, bvh.shape.x);
  ASSERT_EQ(bf.shape.y, bvh.shape.y);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-5f)
          << "Mismatch with noise warp at (" << i << ", " << j << ")";
      EXPECT_NEAR(def(i, j), bvh(i, j), 1e-5f);
      EXPECT_NEAR(no_bvh(i, j), bf(i, j), 1e-5f);
    }
  }
}

TEST(BVH2DPathSDF, ComplexMultiSegmentPathEquivalence)
{
  const int              num_points = 150;
  std::vector<glm::vec3> pts;
  pts.reserve(num_points);

  for (int i = 0; i < num_points; ++i)
  {
    float t = float(i) / float(num_points - 1);
    float x = 0.1f + 0.8f * t;
    float y = 0.5f + 0.35f * std::sin(t * 6.2831853f * 4.f);
    pts.emplace_back(x, y, 0.f);
  }

  Path path(pts);
  path.set_closed(false);

  glm::ivec2 shape(48, 48);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array bf = path_sdf_to_array_brute_force(path, shape, bbox);
  Array bvh = path_sdf_to_array_bvh(path, shape, bbox);
  Array def = path_sdf_to_array(path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-4f)
          << "Mismatch on complex path at (" << i << ", " << j << ")";
      EXPECT_NEAR(def(i, j), bvh(i, j), 1e-4f);
    }
  }
}
