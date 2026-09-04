#include <cstdio>
#include <filesystem>

#include <gtest/gtest.h>

#include "highmap.hpp"

using namespace hmap;

// Regression: the TRI_OPTIMIZED (Delaunay) mesh path wrote texture
// coordinates into an unallocated aiMesh::mTextureCoords[0] -> SIGSEGV.
TEST(ExportAsset, TriOptimizedExportsWithoutCrash)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise(NoiseType::SIMPLEX2, shape, {2.f, 2.f}, 1);
  remap(z);

  std::filesystem::path fname = std::filesystem::temp_directory_path() /
                                "hmap_test_tri_optimized";

  // export_asset appends the format extension itself
  std::filesystem::path out = fname.string() + ".glb";

  bool ok = export_asset(fname.string(),
                         z,
                         MeshType::TRI_OPTIMIZED,
                         AssetExportFormat::GLB2,
                         0.2f,
                         "",
                         "",
                         1e-2f);

  EXPECT_TRUE(ok);
  EXPECT_TRUE(std::filesystem::exists(out));
  std::filesystem::remove(out);
}

TEST(ExportAsset, TriExportsWithoutCrash)
{
  glm::ivec2 shape = {64, 64};
  Array      z = noise(NoiseType::SIMPLEX2, shape, {2.f, 2.f}, 1);
  remap(z);

  std::filesystem::path fname = std::filesystem::temp_directory_path() /
                                "hmap_test_tri";

  // export_asset appends the format extension itself
  std::filesystem::path out = fname.string() + ".glb";

  bool ok = export_asset(fname.string(),
                         z,
                         MeshType::TRI,
                         AssetExportFormat::GLB2,
                         0.2f,
                         "",
                         "",
                         1e-2f);

  EXPECT_TRUE(ok);
  EXPECT_TRUE(std::filesystem::exists(out));
  std::filesystem::remove(out);
}
