#include <vector>

#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);

  glm::ivec2 ij_start = {40, 40};
  glm::ivec2 ij_end = {230, 230};
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};
  float      offset_ratio = 0.3f; // in [0, 1]

  // smallest resolution
  hmap::Path path1 = hmap::find_path_midpoint(z,
                                              ij_start,
                                              ij_end,
                                              bbox,
                                              offset_ratio);

  // limit iterations
  int        max_it = 4;
  hmap::Path path2 = hmap::find_path_midpoint(z,
                                              ij_start,
                                              ij_end,
                                              bbox,
                                              offset_ratio,
                                              max_it);

  // spline interpolation
  hmap::Path path3 = hmap::bspline(path2);

  // --- export path to a png file

  hmap::Array w1 = path1.to_array(shape);
  hmap::Array w2 = path2.to_array(shape);
  hmap::Array w3 = path3.to_array(shape);

  hmap::export_banner_png("ex_find_path_midpoint.png",
                          {z, w1, w2, w3},
                          hmap::Cmap::INFERNO);
}
