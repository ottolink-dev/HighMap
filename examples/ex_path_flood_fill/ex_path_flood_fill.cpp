#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 1;

  hmap::BBox2D bbox = {1.f, 2.f, -0.5f, 0.5f};
  hmap::Path   path = hmap::Path(5,
                               seed,
                               hmap::BBox2D({1.2f, 1.8f, -0.3f, 0.3f}));
  path.reorder_nns();
  path.set_closed(true);

  auto z1 = hmap::Array(shape);
  path.resample_by_spacing(0.1f);
  path = hmap::fractalize(path, 8, seed);
  path.to_array(z1, bbox);

  auto z2 = z1;
  path.to_array(z2, bbox, true);

  hmap::export_banner_png("ex_path_flood_fill.png",
                          {z1, z2},
                          hmap::Cmap::INFERNO);
}
