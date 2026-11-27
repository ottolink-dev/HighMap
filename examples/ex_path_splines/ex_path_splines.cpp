#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 6;

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path path = hmap::Path(10, seed, hmap::BBox2D(0.1f, 0.9f, 0.1f, 0.9f));
  path.reorder_nns();

  path.set_closed(true);

  auto z0 = path.to_array(shape, bbox);
  auto z1 = hmap::bezier(path).to_array(shape, bbox);
  auto z2 = hmap::bezier_round(path).to_array(shape, bbox);
  auto z3 = hmap::bspline(path).to_array(shape, bbox);
  auto z4 = hmap::catmullrom(path).to_array(shape, bbox);
  auto z5 = hmap::decasteljau(path).to_array(shape, bbox);

  hmap::export_banner_png("ex_path_splines.png",
                          {z0, z1, z2, z3, z4, z5},
                          hmap::Cmap::INFERNO);
}
