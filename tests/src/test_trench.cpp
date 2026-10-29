#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(TrenchCarving, EmptyInputsHandledGracefully)
{
  hmap::Array z_empty;
  hmap::Path  path;
  path.emplace_back(0.5f, 0.5f, 0.f);

  hmap::trench(z_empty, path, 0.1f);
  EXPECT_TRUE(z_empty.vector.empty());

  hmap::Array z({32, 32}, 1.f);
  hmap::Path  empty_path;
  hmap::trench(z, empty_path, 0.1f);
  EXPECT_EQ(z.shape.x, 32);
}

TEST(TrenchCarving, SinglePointCarvingModifiesCenter)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z(shape, 1.f);
  hmap::Path  path;
  path.emplace_back(0.5f, 0.5f, 0.f);

  hmap::Array mask;
  hmap::trench(z,
               path,
               0.2f,
               /* enable_width_depth_scaling */ false,
               /* enable_width_distance_scaling */ false,
               /* enable_width_curvature_scaling */ false,
               1.f,
               0.5f,
               2.f,
               hmap::RadialProfile::RP_LINEAR,
               1.f,
               hmap::ElevationLongitudinalProfile::ELP_UNCHANGED,
               0.f,
               0.f,
               0.f,
               0.001f,
               1,
               nullptr,
               &mask);

  // center should be carved to ~0.0
  EXPECT_LT(z(32, 32), 0.5f);
  // corner should remain unchanged
  EXPECT_FLOAT_EQ(z(0, 0), 1.f);
  EXPECT_GT(mask(32, 32), 0.5f);
  EXPECT_FLOAT_EQ(mask(0, 0), 0.f);
}

TEST(TrenchCarving, PolylinePathCarvesSmoothTrench)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z(shape, 1.f);

  hmap::Path path;
  path.emplace_back(0.1f, 0.5f, 0.2f);
  path.emplace_back(0.5f, 0.5f, 0.2f);
  path.emplace_back(0.9f, 0.5f, 0.2f);

  hmap::Array mask;
  hmap::trench(z,
               path,
               0.1f,
               /* enable_width_depth_scaling */ false,
               /* enable_width_distance_scaling */ false,
               /* enable_width_curvature_scaling */ false,
               1.f,
               0.5f,
               2.f,
               hmap::RadialProfile::RP_LINEAR,
               1.f,
               hmap::ElevationLongitudinalProfile::ELP_UNCHANGED,
               0.f,
               0.f,
               0.f,
               0.001f,
               1,
               nullptr,
               &mask);

  // points along the horizontal centerline should be carved to ~0.2
  EXPECT_NEAR(z(64, 64), 0.2f, 0.05f);
  EXPECT_NEAR(z(32, 64), 0.2f, 0.05f);
  EXPECT_NEAR(z(96, 64), 0.2f, 0.05f);

  // top and bottom edges should be unaffected
  EXPECT_FLOAT_EQ(z(64, 0), 1.f);
  EXPECT_FLOAT_EQ(z(64, 127), 1.f);
}

TEST(TrenchCarving, CurvatureScalingAndNoise)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z(shape, 1.f);

  hmap::Path path;
  path.emplace_back(0.1f, 0.2f, 0.1f);
  path.emplace_back(0.5f, 0.8f, 0.1f);
  path.emplace_back(0.9f, 0.2f, 0.1f);

  hmap::Array noise(shape, 0.2f);
  hmap::Array mask;

  hmap::trench(z,
               path,
               0.1f,
               /* enable_width_depth_scaling */ true,
               /* enable_width_distance_scaling */ true,
               /* enable_width_curvature_scaling */ true,
               0.5f,
               0.5f,
               2.f,
               hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
               2.f,
               hmap::ElevationLongitudinalProfile::ELP_DECREASING,
               -0.1f,
               0.1f,
               0.1f,
               0.001f,
               1,
               &noise,
               &mask);

  // ensure all values are finite
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_FALSE(std::isnan(z(i, j)));
      EXPECT_FALSE(std::isinf(z(i, j)));
    }
  }
}

TEST(TrenchCarving, VirtualArrayVsArrayComparisonSingleTile)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z_arr(shape, 1.f);

  hmap::Path path;
  path.emplace_back(0.1f, 0.2f, 0.2f);
  path.emplace_back(0.5f, 0.5f, 0.1f);
  path.emplace_back(0.9f, 0.8f, 0.05f);

  float width = 0.08f;

  // Carve on monolithic array
  hmap::trench(z_arr,
               path,
               width,
               /* enable_width_depth_scaling */ false,
               /* enable_width_distance_scaling */ false,
               /* enable_width_curvature_scaling */ false,
               1.f,
               0.5f,
               2.f,
               hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
               2.f,
               hmap::ElevationLongitudinalProfile::ELP_DECREASING,
               0.f,
               0.f,
               0.f,
               0.001f,
               4,
               nullptr,
               nullptr,
               /* tile_size */ 32,
               {0.f, 1.f, 0.f, 1.f});

  // Carve on 1x1 VirtualArray
  hmap::ComputeMode  cm{.mode = hmap::ForEachMode::VA_SEQUENTIAL};
  hmap::VirtualArray z_va(shape, shape, 0, hmap::StorageMode::VA_RAM);
  z_va.fill(1.f, cm);

  hmap::for_each_tile(
      z_va,
      [&](hmap::Array &tile, const hmap::TileRegion &region)
      {
        hmap::trench(tile,
                     path,
                     width,
                     /* enable_width_depth_scaling */ false,
                     /* enable_width_distance_scaling */ false,
                     /* enable_width_curvature_scaling */ false,
                     1.f,
                     0.5f,
                     2.f,
                     hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
                     2.f,
                     hmap::ElevationLongitudinalProfile::ELP_DECREASING,
                     0.f,
                     0.f,
                     0.f,
                     0.001f,
                     4,
                     nullptr,
                     nullptr,
                     /* tile_size */ 32,
                     region.bbox);
      },
      cm);

  hmap::Array z_va_arr = z_va.to_array(cm);

  float max_diff = 0.f;
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      float diff = std::abs(z_arr(i, j) - z_va_arr(i, j));
      max_diff = std::max(max_diff, diff);
    }
  }

  std::cout << "[DEBUG] SingleTile max_diff = " << max_diff << std::endl;
  EXPECT_NEAR(max_diff, 0.f, 1e-5f);
}

TEST(TrenchCarving, VirtualArrayVsArrayComparisonMultiTile)
{
  glm::ivec2  shape = {128, 128};
  hmap::Array z_arr(shape, 1.f);

  hmap::Path path;
  path.emplace_back(0.1f, 0.2f, 0.2f);
  path.emplace_back(0.5f, 0.5f, 0.1f);
  path.emplace_back(0.9f, 0.8f, 0.05f);

  float width = 0.08f;

  // Carve on monolithic array
  hmap::trench(z_arr,
               path,
               width,
               /* enable_width_depth_scaling */ false,
               /* enable_width_distance_scaling */ false,
               /* enable_width_curvature_scaling */ false,
               1.f,
               0.5f,
               2.f,
               hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
               2.f,
               hmap::ElevationLongitudinalProfile::ELP_DECREASING,
               0.f,
               0.f,
               0.f,
               0.001f,
               4,
               nullptr,
               nullptr,
               /* tile_size */ 32,
               {0.f, 1.f, 0.f, 1.f});

  // Carve on 2x2 VirtualArray with halo 16
  hmap::ComputeMode  cm{.mode = hmap::ForEachMode::VA_SEQUENTIAL};
  hmap::VirtualArray z_va(shape, {64, 64}, 16, hmap::StorageMode::VA_RAM);
  z_va.fill(1.f, cm);

  hmap::for_each_tile(
      z_va,
      [&](hmap::Array &tile, const hmap::TileRegion &region)
      {
        hmap::trench(tile,
                     path,
                     width,
                     /* enable_width_depth_scaling */ false,
                     /* enable_width_distance_scaling */ false,
                     /* enable_width_curvature_scaling */ false,
                     1.f,
                     0.5f,
                     2.f,
                     hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
                     2.f,
                     hmap::ElevationLongitudinalProfile::ELP_DECREASING,
                     0.f,
                     0.f,
                     0.f,
                     0.001f,
                     4,
                     nullptr,
                     nullptr,
                     /* tile_size */ 32,
                     region.bbox);
      },
      cm);

  hmap::Array z_va_arr = z_va.to_array(cm);

  float max_diff = 0.f;
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      float diff = std::abs(z_arr(i, j) - z_va_arr(i, j));
      max_diff = std::max(max_diff, diff);
    }
  }

  std::cout << "[DEBUG] MultiTile max_diff = " << max_diff << std::endl;
  EXPECT_NEAR(max_diff, 0.f, 1e-5f);
}

TEST(TrenchCarving, PathReorderingByLengthDescending)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z1(shape, 1.f);
  hmap::Array z2(shape, 1.f);

  // trunk path: long path crossing the domain from (0.1, 0.5) to (0.9, 0.5)
  hmap::Path trunk;
  trunk.emplace_back(0.1f, 0.5f, 0.f);
  trunk.emplace_back(0.5f, 0.5f, 0.f);
  trunk.emplace_back(0.9f, 0.5f, 0.f);

  // tributary: shorter path ending on trunk at (0.5, 0.5)
  hmap::Path tributary;
  tributary.emplace_back(0.5f, 0.2f, 0.f);
  tributary.emplace_back(0.5f, 0.5f, 0.f);

  // pass order: [tributary, trunk]
  std::vector<hmap::Path> paths_short_first = {tributary, trunk};
  // pass order: [trunk, tributary]
  std::vector<hmap::Path> paths_long_first = {trunk, tributary};

  hmap::trench(z1,
               paths_short_first,
               0.1f,
               /* enable_width_depth_scaling */ false,
               /* enable_width_distance_scaling */ true);

  hmap::trench(z2,
               paths_long_first,
               0.1f,
               /* enable_width_depth_scaling */ false,
               /* enable_width_distance_scaling */ true);

  // regardless of input order, trench should produce identical results because
  // paths are reordered by cumulative length descending
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(z1(i, j), z2(i, j), 1e-6f);
    }
  }

  // tributary terminal width should be non-zero at junction (0.5, 0.5)
  EXPECT_LT(z1(32, 32), 0.5f);
}

TEST(TrenchCarving, GraphFlowNetworkWithOutletDistanceScaling)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z(shape, 1.0f);

  // create Y junction graph with endpoints placed on grid coordinates
  // (x = i / 64, y = j / 64)
  int i_root0 = 13;
  int j_root0 = 51;
  int i_root1 = 51;
  int j_root1 = 51;
  int i_conf = 32;
  int j_conf = 32;
  int i_mouth = 32;
  int j_mouth = 13;

  hmap::Graph g;
  g.push_back(hmap::Point(float(i_root0) / 64.f,
                          float(j_root0) / 64.f,
                          1.0f)); // 0: branch 1 head
  g.push_back(hmap::Point(float(i_root1) / 64.f,
                          float(j_root1) / 64.f,
                          1.0f)); // 1: branch 2 head
  g.push_back(hmap::Point(float(i_conf) / 64.f,
                          float(j_conf) / 64.f,
                          0.5f)); // 2: confluence
  g.push_back(hmap::Point(float(i_mouth) / 64.f,
                          float(j_mouth) / 64.f,
                          0.1f)); // 3: trunk mouth

  g.add_edge(0, 2);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  // test compute_relative_widths
  std::vector<float> rel_w = g.compute_relative_widths(1.0f);
  ASSERT_EQ(rel_w.size(), 4u);
  EXPECT_FLOAT_EQ(rel_w[0], 0.f);
  EXPECT_FLOAT_EQ(rel_w[1], 0.f);
  EXPECT_GT(rel_w[2], 0.f);
  EXPECT_LT(rel_w[2], 1.f);
  EXPECT_FLOAT_EQ(rel_w[3], 1.f);

  float target_width = 0.08f;

  // carve with explicit relative_widths
  hmap::trench(z,
               g,
               rel_w,
               target_width,
               /* enable_width_depth_scaling */ false);

  // mouth (outlet) is carved
  EXPECT_LT(z(i_mouth, j_mouth), 1.0f);

  // mouth at outlet should have full width incision:
  // at lateral offset 0.035 (~2 cells) from mouth center (32, 13), it should be
  // carved
  int i_mouth_offset = i_mouth + static_cast<int>(std::round(0.035f * 64.f));
  EXPECT_LT(z(i_mouth_offset, j_mouth), 1.0f);
}

TEST(CarveRiver, RiverNetworkIncisionAndMonotonicity)
{
  glm::ivec2  shape = {64, 64};
  hmap::Array z(shape, 1.0f);
  hmap::Array z_orig = z;

  // Y-confluence network
  hmap::Graph g;
  g.push_back(hmap::Point(0.2f, 0.8f, 1.0f)); // 0: headwater 1
  g.push_back(hmap::Point(0.8f, 0.8f, 0.9f)); // 1: headwater 2
  g.push_back(hmap::Point(0.5f, 0.5f, 0.6f)); // 2: confluence
  g.push_back(hmap::Point(0.5f, 0.1f, 0.2f)); // 3: outlet

  g.add_edge(0, 2);
  g.add_edge(1, 2);
  g.add_edge(2, 3);

  float outlet_width = 0.10f;
  float headwater_width = 0.02f;
  float depth = 0.15f;

  hmap::Array mask;
  hmap::carve_river(z,
                    g,
                    outlet_width,
                    headwater_width,
                    /* width_exponent */ 0.5f,
                    depth,
                    /* bank_width_ratio */ 1.0f,
                    hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
                    2.0f,
                    /* p_noise_r */ nullptr,
                    &mask);

  // 1. Strict incision: terrain should never be raised
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_LE(z(i, j), z_orig(i, j) + 1e-6f);
      EXPECT_FALSE(std::isnan(z(i, j)));
      EXPECT_FALSE(std::isinf(z(i, j)));
      EXPECT_GE(mask(i, j), 0.f);
      EXPECT_LE(mask(i, j), 1.f);
    }
  }

  // 2. Channel floor incision along stream path
  // At the outlet node (0.5, 0.1) -> pixel (32, 6)
  int i_out = static_cast<int>(std::round(0.5f * (shape.x - 1)));
  int j_out = static_cast<int>(std::round(0.1f * (shape.y - 1)));
  EXPECT_LT(z(i_out, j_out), 1.0f);
  EXPECT_NEAR(z(i_out, j_out), 0.2f - depth, 0.05f);

  // 3. Tributary confluence is carved
  int i_conf = static_cast<int>(std::round(0.5f * (shape.x - 1)));
  int j_conf = static_cast<int>(std::round(0.5f * (shape.y - 1)));
  EXPECT_LT(z(i_conf, j_conf), 1.0f);

  // 4. Outlet has wider incision than headwater
  // Measure carved width along horizontal slices
  int carved_cells_headwater = 0;
  int j_head = static_cast<int>(std::round(0.8f * (shape.y - 1)));
  for (int i = 0; i < shape.x; ++i)
  {
    if (mask(i, j_head) > 0.05f) carved_cells_headwater++;
  }

  int carved_cells_outlet = 0;
  for (int i = 0; i < shape.x; ++i)
  {
    if (mask(i, j_out) > 0.05f) carved_cells_outlet++;
  }

  EXPECT_GT(carved_cells_outlet, carved_cells_headwater);
}
