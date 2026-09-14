#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {2.f, 2.f};
  int        seed = 0;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);

  // Dijkstra
  hmap::Path path = hmap::find_cut_path_dijkstra(
      z,
      hmap::DomainBoundary::BOUNDARY_LEFT,
      hmap::DomainBoundary::BOUNDARY_RIGHT);
  hmap::Array zp = path.to_array(shape);

  // greedy midpoint
  hmap::Path path2 = hmap::find_cut_path_midpoint(
      z,
      hmap::DomainBoundary::BOUNDARY_LEFT,
      hmap::DomainBoundary::BOUNDARY_RIGHT,
      seed);
  hmap::Array zp2 = path2.to_array(shape);

  // multiscale
  hmap::Path path3 = hmap::find_cut_path_multiscale(
      z,
      hmap::DomainBoundary::BOUNDARY_LEFT,
      hmap::DomainBoundary::BOUNDARY_RIGHT);
  hmap::Array zp3 = path3.to_array(shape);

  hmap::export_banner_png("ex_find_cut_path.png",
                          {z, zp, zp2, zp3},
                          hmap::Cmap::INFERNO);
}
