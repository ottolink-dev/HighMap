#include <vector>

#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  z.to_png("ex_find_path_dijkstra0.png", hmap::Cmap::TERRAIN, true);

  glm::ivec2 ij_start = {40, 40};
  glm::ivec2 ij_end = {230, 230};
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  // default parameters
  hmap::Path  path1 = hmap::find_path_dijkstra(z, ij_start, ij_end, bbox);
  hmap::Array w1 = path1.to_array(shape);

  // set "elevation_ratio" to 1.f to find the path with the lowest
  // cumulative elevation
  hmap::Path  path2 = hmap::find_path_dijkstra(z, ij_start, ij_end, bbox, 1.f);
  hmap::Array w2 = path2.to_array(shape);

  hmap::export_banner_png("ex_find_path_dijkstra.png",
                          {z, w1, w2},
                          hmap::Cmap::TERRAIN);
}
