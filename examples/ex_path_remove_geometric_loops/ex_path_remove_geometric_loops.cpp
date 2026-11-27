#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 0;

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path   path = hmap::Path(40, seed, bbox);
  path.set_closed(false);

  auto z1 = hmap::Array(shape);
  path.to_array(z1, bbox);

  auto       z2 = hmap::Array(shape);
  hmap::Path path_noloop = hmap::remove_geometric_loops(path);
  path_noloop.to_array(z2, bbox);

  hmap::export_banner_png("ex_path_remove_geometry_loops.png",
                          {z1, z2},
                          hmap::Cmap::INFERNO);
}
