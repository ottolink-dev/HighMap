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
