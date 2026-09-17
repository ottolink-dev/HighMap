#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap.hpp"
#include "highmap/hydrology/drainage_basin.hpp"

#include <gtest/gtest.h>

TEST(DrainageBasin, BasicConstructionAndOutlets)
{
  std::vector<glm::vec3> pts = {{0.f, 0.f, 0.f},
                                {1.f, 0.f, 0.f},
                                {1.f, 1.f, 0.f},
                                {0.f, 1.f, 0.f},
                                {0.5f, 0.5f, 1.f}};

  hmap::DrainageBasin db(pts);
  EXPECT_EQ(db.size(), 5);

  db.update_stream_tree();
  const auto &outlets = db.get_outlets();
  EXPECT_FALSE(outlets.empty());
}

TEST(DrainageBasin, StrahlerOrderTree)
{
  std::vector<glm::vec3> pts = {{0.f, 0.f, 0.f},
                                {1.f, 0.f, 1.f},
                                {1.f, 1.f, 2.f},
                                {0.f, 1.f, 1.f},
                                {0.5f, 0.5f, 1.5f}};

  hmap::DrainageBasin db(pts);
  db.set_outlets({0});
  db.update_stream_tree();

  std::vector<size_t> order = db.compute_strahler_order();
  EXPECT_EQ(order.size(), pts.size());

  for (size_t o : order)
    EXPECT_GE(o, 1);
}

TEST(DrainageBasin, RidgeNodesAndRootsInitialized)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   {2.f, 2.f},
                                   42);
  hmap::remap(z0);

  glm::vec4   bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Cloud cloud = hmap::random_cloud(200,
                                         42,
                                         hmap::PointSamplingMethod::RND_LHS,
                                         bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z0, bbox);

  hmap::DrainageBasin db(cloud.to_vec3());
  db.update_stream_tree();

  std::vector<bool> ridges = db.compute_is_ridge_node();
  EXPECT_EQ(ridges.size(), db.size());
}

TEST(DrainageBasin, ResponseTimesAndElevationUpdate)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   {2.f, 2.f},
                                   123);
  hmap::remap(z0);

  glm::vec4   bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Cloud cloud = hmap::random_cloud(300,
                                         123,
                                         hmap::PointSamplingMethod::RND_LHS,
                                         bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z0, bbox);

  hmap::DrainageBasin db(cloud.to_vec3());
  db.set_outlets(hmap::find_border_sinks(db.get_mesh()));
  db.update_stream_tree();

  std::vector<float> area = db.compute_vertex_areas();
  std::vector<float> acc(db.size(), 0.f);
  db.accumulate_area_by_outlet(area, acc);

  std::vector<float> erodibility(db.size(), 1.f);
  std::vector<float> response_times = db.compute_response_times(acc,
                                                                erodibility,
                                                                0.8f);
  EXPECT_EQ(response_times.size(), db.size());

  std::vector<float> max_slope(db.size(), 2.f);
  float diff = db.update_elevations(response_times, 1.f, max_slope);
  EXPECT_GE(diff, 0.f);
}
