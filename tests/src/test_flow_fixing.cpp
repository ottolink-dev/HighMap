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
