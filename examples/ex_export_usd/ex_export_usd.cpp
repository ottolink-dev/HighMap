#include <iostream>
#include <vector>

#include "highmap.hpp"

int main(void)
{
  // --- Terrain heightmap

  glm::ivec2 shape = {128, 128};
  glm::vec2  kv = {2.f, 2.f};
  int        seed = 42;

  hmap::Array terrain = hmap::noise(hmap::NoiseType::SIMPLEX2, shape, kv, seed);
  hmap::remap(terrain);

  // --- Clouds (scattered 3D point features)

  hmap::Cloud              trees(50, seed + 1);
  hmap::Cloud              rocks(30, seed + 2);
  std::vector<hmap::Cloud> clouds = {trees, rocks};

  // --- Paths (polylines / roads / rivers)

  std::vector<hmap::Point> road_points = {hmap::Point(0.1f, 0.1f, 0.2f),
                                          hmap::Point(0.3f, 0.4f, 0.35f),
                                          hmap::Point(0.6f, 0.5f, 0.45f),
                                          hmap::Point(0.9f, 0.85f, 0.25f)};
  hmap::Path               road(road_points);

  std::vector<hmap::Point> river_points = {hmap::Point(0.05f, 0.85f, 0.15f),
                                           hmap::Point(0.45f, 0.55f, 0.22f),
                                           hmap::Point(0.85f, 0.15f, 0.1f)};
  hmap::Path               river(river_points);

  std::vector<hmap::Path> paths = {road, river};

  // --- Export composite scenes to USD formats

  std::cout << "Exporting to scene.usda (ASCII USD)...\n";
  hmap::export_usd("scene.usda",
                   terrain,
                   clouds,
                   paths,
                   hmap::MeshType::TRI,
                   0.25f);

  std::cout << "Exporting to scene.usdc (Binary Crate USD)...\n";
  hmap::export_usd("scene.usdc",
                   terrain,
                   clouds,
                   paths,
                   hmap::MeshType::TRI,
                   0.25f);

  std::cout << "Exporting to scene_opt.usda (Delaunay optimized terrain)...\n";
  hmap::export_usd("scene_opt.usda",
                   terrain,
                   clouds,
                   paths,
                   hmap::MeshType::TRI_OPTIMIZED,
                   0.25f,
                   "",
                   "",
                   1e-2f);

  std::cout << "USD export completed successfully.\n";
  return 0;
}
