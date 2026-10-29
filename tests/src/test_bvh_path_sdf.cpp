#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/graph.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/sdf.hpp"

#include <gtest/gtest.h>

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

TEST(BVH2DPathSDF, BezierMatchesBruteForceAndBVH)
{
  Path path(std::vector<Point>{{0.1f, 0.2f, 0.f},
                               {0.3f, 0.8f, 0.f},
                               {0.5f, 0.2f, 0.f},
                               {0.7f, 0.9f, 0.f},
                               {0.9f, 0.3f, 0.f}});
  path.set_closed(false);

  glm::ivec2 shape(32, 32);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array bf = sdf_2d_polyline_bezier_brute_force(path, shape, bbox);
  Array bvh = sdf_2d_polyline_bezier_bvh(path, shape, bbox);
  Array def = sdf_2d_polyline_bezier(path, shape, bbox);
  Array no_bvh = sdf_2d_polyline_bezier(path,
                                        shape,
                                        bbox,
                                        nullptr,
                                        nullptr,
                                        false);

  ASSERT_EQ(bf.shape.x, bvh.shape.x);
  ASSERT_EQ(bf.shape.y, bvh.shape.y);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(bf(i, j)));
      EXPECT_FALSE(std::isnan(bvh(i, j)));
      EXPECT_FALSE(std::isinf(bf(i, j)));
      EXPECT_FALSE(std::isinf(bvh(i, j)));

      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-4f)
          << "Bezier mismatch at pixel (" << i << ", " << j << ")";
      EXPECT_NEAR(def(i, j), bvh(i, j), 1e-4f);
      EXPECT_NEAR(no_bvh(i, j), bf(i, j), 1e-4f);
    }
  }
}

TEST(BVH2DPathSDF, BezierColinearPointsNoNaN)
{
  Path path(std::vector<Point>{{0.1f, 0.1f, 0.f},
                               {0.5f, 0.5f, 0.f},
                               {0.9f, 0.9f, 0.f}});
  path.set_closed(false);

  glm::ivec2 shape(16, 16);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array res_bvh = sdf_2d_polyline_bezier_bvh(path, shape, bbox);
  Array res_bf = sdf_2d_polyline_bezier_brute_force(path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(res_bvh(i, j)));
      EXPECT_FALSE(std::isinf(res_bvh(i, j)));
      EXPECT_FALSE(std::isnan(res_bf(i, j)));
      EXPECT_FALSE(std::isinf(res_bf(i, j)));
      EXPECT_NEAR(res_bvh(i, j), res_bf(i, j), 1e-5f);
    }
  }
}

TEST(BVH2DPathSDF, BezierClosedPathSigned)
{
  Path path(std::vector<Point>{{0.25f, 0.25f, 0.f},
                               {0.75f, 0.25f, 0.f},
                               {0.75f, 0.75f, 0.f},
                               {0.25f, 0.75f, 0.f}});
  path.set_closed(true);

  glm::ivec2 shape(24, 24);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array bf = sdf_2d_polyline_bezier_brute_force(path, shape, bbox);
  Array bvh = sdf_2d_polyline_bezier_bvh(path, shape, bbox);
  Array def = sdf_2d_polyline_bezier(path, shape, bbox);

  ASSERT_EQ(bf.shape.x, bvh.shape.x);
  ASSERT_EQ(bf.shape.y, bvh.shape.y);

  // Closed path must have negative values inside and positive values outside
  EXPECT_LT(bvh.min(), 0.f);
  EXPECT_GT(bvh.max(), 0.f);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(bvh(i, j)));
      EXPECT_FALSE(std::isinf(bvh(i, j)));
      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-4f)
          << "Closed Bezier SDF mismatch at (" << i << ", " << j << ")";
      EXPECT_NEAR(def(i, j), bvh(i, j), 1e-4f);
    }
  }
}

TEST(BVH2DPathSDF, BezierEvenPointCountHandled)
{
  Path path(std::vector<Point>{{0.1f, 0.1f, 0.f},
                               {0.3f, 0.7f, 0.f},
                               {0.6f, 0.8f, 0.f},
                               {0.9f, 0.2f, 0.f}});
  path.set_closed(false);

  glm::ivec2 shape(20, 20);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array res_bvh = sdf_2d_polyline_bezier_bvh(path, shape, bbox);
  Array res_bf = sdf_2d_polyline_bezier_brute_force(path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(res_bvh(i, j)));
      EXPECT_FALSE(std::isinf(res_bvh(i, j)));
      EXPECT_NEAR(res_bvh(i, j), res_bf(i, j), 1e-4f);
    }
  }
}

TEST(BVH2DPathSDF, GpuBezierMatchesCpu)
{
  // Test open path
  Path open_path(std::vector<Point>{{0.1f, 0.2f, 0.f},
                                    {0.4f, 0.8f, 0.f},
                                    {0.7f, 0.3f, 0.f},
                                    {0.9f, 0.6f, 0.f}});
  open_path.set_closed(false);

  glm::ivec2 shape(32, 32);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array cpu_open = sdf_2d_polyline_bezier(open_path, shape, bbox);
  Array gpu_open = gpu::sdf_2d_polyline_bezier(open_path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(gpu_open(i, j)));
      EXPECT_FALSE(std::isinf(gpu_open(i, j)));
      EXPECT_NEAR(gpu_open(i, j), cpu_open(i, j), 1e-3f)
          << "GPU open Bezier mismatch at (" << i << ", " << j << ")";
    }
  }

  // Test colinear path (division-by-zero check on GPU)
  Path colinear_path(std::vector<Point>{{0.1f, 0.1f, 0.f},
                                        {0.5f, 0.5f, 0.f},
                                        {0.9f, 0.9f, 0.f}});
  colinear_path.set_closed(false);

  Array cpu_colinear = sdf_2d_polyline_bezier(colinear_path, shape, bbox);
  Array gpu_colinear = gpu::sdf_2d_polyline_bezier(colinear_path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(gpu_colinear(i, j)));
      EXPECT_FALSE(std::isinf(gpu_colinear(i, j)));
      EXPECT_NEAR(gpu_colinear(i, j), cpu_colinear(i, j), 1e-3f)
          << "GPU colinear Bezier mismatch at (" << i << ", " << j << ")";
    }
  }

  // Test closed path (sign / ray-crossing parity check on GPU)
  Path closed_path(std::vector<Point>{{0.25f, 0.25f, 0.f},
                                      {0.75f, 0.25f, 0.f},
                                      {0.75f, 0.75f, 0.f},
                                      {0.25f, 0.75f, 0.f}});
  closed_path.set_closed(true);

  Array cpu_closed = sdf_2d_polyline_bezier(closed_path, shape, bbox);
  Array gpu_closed = gpu::sdf_2d_polyline_bezier(closed_path, shape, bbox);

  EXPECT_LT(gpu_closed.min(), 0.f);
  EXPECT_GT(gpu_closed.max(), 0.f);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(gpu_closed(i, j)));
      EXPECT_FALSE(std::isinf(gpu_closed(i, j)));
      EXPECT_NEAR(gpu_closed(i, j), cpu_closed(i, j), 1e-3f)
          << "GPU closed Bezier mismatch at (" << i << ", " << j << ")";
    }
  }
}

TEST(BVH2DPathSDF, GpuPathSdfToArrayMatchesCpu)
{
  // Test open path
  Path open_path(std::vector<Point>{{0.2f, 0.3f, 0.f},
                                    {0.5f, 0.8f, 0.f},
                                    {0.8f, 0.4f, 0.f},
                                    {0.9f, 0.7f, 0.f}});
  open_path.set_closed(false);

  glm::ivec2 shape(32, 32);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array cpu_open = path_sdf_to_array(open_path, shape, bbox);
  Array gpu_open = gpu::path_sdf_to_array(open_path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(gpu_open(i, j)));
      EXPECT_FALSE(std::isinf(gpu_open(i, j)));
      EXPECT_NEAR(gpu_open(i, j), cpu_open(i, j), 1e-3f)
          << "GPU open path_sdf_to_array mismatch at (" << i << ", " << j
          << ")";
    }
  }

  // Test closed path with interior negative signs
  Path closed_path(std::vector<Point>{{0.25f, 0.25f, 0.f},
                                      {0.75f, 0.25f, 0.f},
                                      {0.75f, 0.75f, 0.f},
                                      {0.25f, 0.75f, 0.f}});
  closed_path.set_closed(true);

  Array cpu_closed = path_sdf_to_array(closed_path, shape, bbox);
  Array gpu_closed = gpu::path_sdf_to_array(closed_path, shape, bbox);

  EXPECT_LT(gpu_closed.min(), 0.f);
  EXPECT_GT(gpu_closed.max(), 0.f);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(gpu_closed(i, j)));
      EXPECT_FALSE(std::isinf(gpu_closed(i, j)));
      EXPECT_NEAR(gpu_closed(i, j), cpu_closed(i, j), 1e-3f)
          << "GPU closed path_sdf_to_array mismatch at (" << i << ", " << j
          << ")";
    }
  }

  // Test closed path with noise domain warping
  Array dx(shape, 0.05f);
  Array dy(shape, -0.03f);

  Array cpu_warped = path_sdf_to_array(closed_path, shape, bbox, &dx, &dy);
  Array gpu_warped = gpu::path_sdf_to_array(closed_path, shape, bbox, &dx, &dy);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(gpu_warped(i, j)));
      EXPECT_FALSE(std::isinf(gpu_warped(i, j)));
      EXPECT_NEAR(gpu_warped(i, j), cpu_warped(i, j), 1e-3f)
          << "GPU warped path_sdf_to_array mismatch at (" << i << ", " << j
          << ")";
    }
  }
}

TEST(BVH2DPathSDF, PathSdfBezierToArrayMatchesSdf2DPolylineBezier)
{
  Path path(std::vector<Point>{{0.1f, 0.2f, 0.f},
                               {0.3f, 0.8f, 0.f},
                               {0.6f, 0.1f, 0.f},
                               {0.9f, 0.7f, 0.f},
                               {0.5f, 0.5f, 0.f}});

  glm::ivec2 shape(32, 32);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array bf = path_sdf_bezier_to_array_brute_force(path, shape, bbox);
  Array bvh = path_sdf_bezier_to_array_bvh(path, shape, bbox);
  Array def = path_sdf_bezier_to_array(path, shape, bbox);
  Array no_bvh = path_sdf_bezier_to_array(path,
                                          shape,
                                          bbox,
                                          nullptr,
                                          nullptr,
                                          false);
  Array ref = sdf_2d_polyline_bezier(path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-4f);
      EXPECT_FLOAT_EQ(def(i, j), bvh(i, j));
      EXPECT_FLOAT_EQ(no_bvh(i, j), bf(i, j));
      EXPECT_FLOAT_EQ(def(i, j), ref(i, j));
    }
  }

  // GPU test
  Array gpu_res = gpu::path_sdf_bezier_to_array(path, shape, bbox);
  Array gpu_ref = gpu::sdf_2d_polyline_bezier(path, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FLOAT_EQ(gpu_res(i, j), gpu_ref(i, j));
      EXPECT_NEAR(gpu_res(i, j), def(i, j), 1e-3f);
    }
  }
}

TEST(BVH2DGraphSDF, MatchesBruteForceAndBVH)
{
  Graph graph;

  graph.push_back(Point(0.2f, 0.2f, 0.f));
  graph.push_back(Point(0.5f, 0.5f, 0.f));
  graph.push_back(Point(0.8f, 0.2f, 0.f));
  graph.push_back(Point(0.5f, 0.8f, 0.f));

  graph.add_edge(0, 1);
  graph.add_edge(2, 1);
  graph.add_edge(1, 3);

  glm::ivec2 shape(32, 32);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array bf = sdf_2d_graph_brute_force(graph, shape, bbox);
  Array bvh = sdf_2d_graph_bvh(graph, shape, bbox);
  Array def = sdf_2d_graph(graph, shape, bbox);
  Array to_arr = graph.to_array_sdf(shape, bbox);

  ASSERT_EQ(bf.shape.x, bvh.shape.x);
  ASSERT_EQ(bf.shape.y, bvh.shape.y);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(bvh(i, j), bf(i, j), 1e-5f)
          << "Mismatch at (" << i << ", " << j << ")";
      EXPECT_FLOAT_EQ(def(i, j), bvh(i, j));
      EXPECT_FLOAT_EQ(to_arr(i, j), bvh(i, j));
    }
  }
}

TEST(BVH2DGraphSDF, EmptyGraphYieldsLargeDistance)
{
  Graph      graph;
  glm::ivec2 shape(16, 16);
  glm::vec4  bbox(0.f, 1.f, 0.f, 1.f);

  Array res = sdf_2d_graph(graph, shape, bbox);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_GE(res(i, j), 1e5f);
    }
  }
}
