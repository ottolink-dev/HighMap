#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 2;
  int        npoints = 10;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path   path = hmap::Path(npoints,
                               seed,
                               hmap::BBox2D(0.1f, 0.9f, 0.1f, 0.9f));
  path.set_closed(false);
  path.reorder_nns();
  path.set_values_from_array(z, bbox);

  // before
  auto z1 = path.to_array(shape, bbox);

  // after
  int edge_divisions = 0;
  path = hmap::dijkstra(path, z, hmap::UvRect::unit(), edge_divisions, 0.9f);
  path.set_values_from_array(z, bbox);

  auto z2 = path.to_array(shape, bbox);

  hmap::export_banner_png("ex_path_dijkstra.png",
                          {z, z1, z2},
                          hmap::Cmap::INFERNO);
}
