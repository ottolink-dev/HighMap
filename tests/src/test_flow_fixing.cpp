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

TEST(FlowFixingMST, VirtualArrayTriangulatedExecution)
{
  glm::ivec2        shape = {128, 128};
  glm::ivec2        tile_shape = {64, 64};
  int               halo = 16;
  hmap::ComputeMode cm;

  hmap::Array z_arr = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {2.f, 2.f},
                                      42,
                                      8,
                                      2.f,
                                      0.5f,
                                      1.f,
                                      nullptr,
                                      nullptr,
                                      nullptr,
                                      {0.f, 1.f, 0.f, 1.f});

  hmap::VirtualArray z_va(shape, tile_shape, halo, hmap::StorageMode::VA_RAM);
  z_va.from_array(z_arr, cm);

  hmap::VirtualArray z_fixed_va = hmap::va::flow_fixing_mst_triangulated(
      cm,
      z_va,
      512,                                      // control_points_count
      42,                                       // seed
      0.001f,                                   // riverbed_talus
      0.95f,                                    // elevation_ratio
      2.f,                                      // distance_exponent
      10.f,                                     // upward_penalization
      1e-4f,                                    // minimum_depth
      4.f,                                      // merging_distance
      hmap::RadialProfile::RP_SMOOTHSTEP_UPPER, // radial_profile
      2.f);

  hmap::Array z_fixed_arr = hmap::flow_fixing_mst_triangulated(
      z_arr,
      512,
      42,
      0.001f,
      0.95f,
      2.f,
      10.f,
      1e-4f,
      4.f,
      hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
      2.f);

  hmap::Array z_va_to_arr = z_fixed_va.to_array(cm);

  float max_diff = 0.f;
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      float d = std::abs(z_fixed_arr(i, j) - z_va_to_arr(i, j));
      max_diff = std::max(max_diff, d);
    }
  }

  EXPECT_NEAR(max_diff, 0.f, 0.1f);
}

TEST(FlowAccumulationDinf, StandaloneVsTiledCommunicator)
{
  glm::ivec2 shape = {64, 64};
  glm::ivec2 tile_shape = {32, 32};
  const int  halo = 4;
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  // Generate test terrain: inclined plane with some perlin noise
  hmap::Array z = hmap::noise(hmap::NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  hmap::remap(z, 0.f, 0.5f);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      // Sloping from (0,0) down to (63,63)
      z(i, j) += 2.0f - float(i + j) / float(shape.x + shape.y);
    }
  }

  // 1. Monolithic standalone reference
  hmap::Array facc_ref = hmap::flow_accumulation_dinf(z, 0.1f);

  EXPECT_EQ(facc_ref.shape.x, shape.x);
  EXPECT_EQ(facc_ref.shape.y, shape.y);
  EXPECT_GT(facc_ref.max(), 1.f);
  EXPECT_GE(facc_ref.min(), 1.f);

  // 2. Distributed execution across tiles using for_each_tile with
  // ICommunicator
  hmap::VirtualArray z_va(shape,
                          bbox,
                          tile_shape,
                          halo,
                          hmap::StorageMode::VA_RAM);
  hmap::VirtualArray facc_va(shape,
                             bbox,
                             tile_shape,
                             halo,
                             hmap::StorageMode::VA_RAM);

  hmap::ComputeMode cm{.mode = hmap::ForEachMode::VA_DISTRIBUTED};
  z_va.from_array(z, cm);

  hmap::for_each_tile(
      {&z_va},
      {&facc_va},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>       p_arrays_out,
          const hmap::TileRegion          &region,
          hmap::ICommunicator             &comm)
      {
        auto [pa_z] = hmap::unpack<1>(p_arrays_in);
        auto [pa_facc] = hmap::unpack<1>(p_arrays_out);

        *pa_facc = hmap::flow_accumulation_dinf(*pa_z, 0.1f, 0, region, comm);
      },
      cm);

  hmap::Array facc_tiled = facc_va.to_array(cm);

  EXPECT_EQ(facc_tiled.shape.x, shape.x);
  EXPECT_EQ(facc_tiled.shape.y, shape.y);
  EXPECT_GE(facc_tiled.min(), 1.f);

  // Flow accumulation should have accumulated significantly in downstream
  // regions
  EXPECT_GT(facc_tiled.max(), 10.f);

  // Mean flow accumulation between monolithic and tiled should correlate well
  EXPECT_NEAR(facc_tiled.mean(), facc_ref.mean(), facc_ref.mean() * 0.35f);
}
