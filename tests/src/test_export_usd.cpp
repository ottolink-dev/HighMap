#include <filesystem>

#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(ExportUsd, BasicExport)
{
  glm::ivec2  shape = {32, 32};
  hmap::Array terrain = hmap::noise(hmap::NoiseType::SIMPLEX2,
                                    shape,
                                    {2.f, 2.f},
                                    1);
  hmap::remap(terrain);

  // Create clouds
  hmap::Cloud              cloud1(10, 42);
  std::vector<hmap::Cloud> clouds = {cloud1};

  // Create path
  std::vector<hmap::Point> pts = {hmap::Point(0.1f, 0.1f, 0.5f),
                                  hmap::Point(0.5f, 0.5f, 0.8f),
                                  hmap::Point(0.9f, 0.9f, 0.2f)};
  hmap::Path               path(pts);
  std::vector<hmap::Path>  paths = {path};

  std::filesystem::path temp_usda = std::filesystem::temp_directory_path() /
                                    "test_highmap_scene.usda";
  std::filesystem::path temp_usdc = std::filesystem::temp_directory_path() /
                                    "test_highmap_scene.usdc";

  // Export USDA (ASCII)
  bool ok_usda = hmap::export_usd(temp_usda.string(), terrain, clouds, paths);
  EXPECT_TRUE(ok_usda);
  EXPECT_TRUE(std::filesystem::exists(temp_usda));
  EXPECT_GT(std::filesystem::file_size(temp_usda), 0u);

  // Export USDC (Binary Crate)
  bool ok_usdc = hmap::export_usd(temp_usdc.string(), terrain, clouds, paths);
  EXPECT_TRUE(ok_usdc);
  EXPECT_TRUE(std::filesystem::exists(temp_usdc));
  EXPECT_GT(std::filesystem::file_size(temp_usdc), 0u);

  // Optimized triangulation export
  std::filesystem::path temp_opt = std::filesystem::temp_directory_path() /
                                   "test_highmap_opt.usda";
  bool ok_opt = hmap::export_usd(temp_opt.string(),
                                 terrain,
                                 clouds,
                                 paths,
                                 hmap::MeshType::TRI_OPTIMIZED);
  EXPECT_TRUE(ok_opt);
  EXPECT_TRUE(std::filesystem::exists(temp_opt));

  // Clean up
  std::filesystem::remove(temp_usda);
  std::filesystem::remove(temp_usdc);
  std::filesystem::remove(temp_opt);
}
