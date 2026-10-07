#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(FlowFixingMST, DijkstraModePreservesShapeAndNonEmpty)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);
  hmap::remap(z);

  float       riverbed_talus = 0.01f / shape.x;
  hmap::Array z_fixed = hmap::flow_fixing_mst(
      z,
      riverbed_talus,
      0.95f,
      2.f,
      10.f,
      0.5f,
      4,
      1e-4f,
      true,
      4.f,
      hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
      2.f,
      nullptr,
      false);

  EXPECT_EQ(z_fixed.shape.x, shape.x);
  EXPECT_EQ(z_fixed.shape.y, shape.y);
}

TEST(FlowFixingMST, MidpointModePreservesShapeAndNonEmpty)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);
  hmap::remap(z);

  float       riverbed_talus = 0.01f / shape.x;
  hmap::Array z_fixed = hmap::flow_fixing_mst(
      z,
      riverbed_talus,
      0.95f,
      2.f,
      10.f,
      0.5f,
      4,
      1e-4f,
      true,
      4.f,
      hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
      2.f,
      nullptr,
      true);

  EXPECT_EQ(z_fixed.shape.x, shape.x);
  EXPECT_EQ(z_fixed.shape.y, shape.y);
}

TEST(FlowFixingMST, MidpointReducesSinks)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {3.f, 3.f},
                                  123);
  hmap::remap(z);

  std::vector<glm::ivec2> sinks_before = hmap::find_flow_sinks(z);

  float       riverbed_talus = 0.01f / shape.x;
  hmap::Array z_fixed = hmap::flow_fixing_mst(
      z,
      riverbed_talus,
      0.95f,
      2.f,
      50.f,
      0.5f,
      4,
      1e-4f,
      true,
      4.f,
      hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
      2.f,
      nullptr,
      true);

  std::vector<glm::ivec2> sinks_after = hmap::find_flow_sinks(z_fixed);

  // sink count should be significantly reduced or equal
  EXPECT_LE(sinks_after.size(), sinks_before.size());
}

TEST(FlowFixingMST, TriangulatedPathsAndCarving)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {2.f, 2.f},
                                  42);
  hmap::remap(z);

  // generate paths via triangulated mesh
  std::vector<hmap::Path> paths = hmap::flow_fixing_mst_paths(z,
                                                              512,
                                                              42,
                                                              0.001f);

  // verify paths structure
  for (const auto &path : paths)
  {
    EXPECT_GE(path.size(), 2);
    // verify coordinates are within [0, 1] bounding box
    for (size_t i = 0; i < path.size(); ++i)
    {
      const auto &pt = path[i];
      EXPECT_GE(pt.x, 0.f);
      EXPECT_LE(pt.x, 1.f);
      EXPECT_GE(pt.y, 0.f);
      EXPECT_LE(pt.y, 1.f);
    }
  }

  // test carving via triangulated flow fixing
  hmap::Array z_tri_fixed = hmap::flow_fixing_mst_triangulated(z,
                                                               512,
                                                               42,
                                                               0.001f,
                                                               0.95f,
                                                               2.f,
                                                               10.f,
                                                               1e-4f,
                                                               4.f);

  EXPECT_EQ(z_tri_fixed.shape.x, shape.x);
  EXPECT_EQ(z_tri_fixed.shape.y, shape.y);
}
