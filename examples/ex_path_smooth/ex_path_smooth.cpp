#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 6;
  int        npoints = 10;

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path   path = hmap::Path(npoints,
                               seed,
                               hmap::BBox2D(0.1f, 0.9f, 0.1f, 0.9f));
  path.reorder_nns();
  path = hmap::fractalize(path, 4, seed, 0.2f);
  path.resample_by_spacing(0.05f);

  auto z1 = path.to_array(shape, bbox);

  int        navg = 12;
  hmap::Path path_c = hmap::smooth(path, navg);
  auto       z2 = path_c.to_array(shape, bbox);

  hmap::assert_start_end_points(path, path_c);

  path.set_closed(true);
  path_c = hmap::smooth(path, navg);
  auto z3 = path_c.to_array(shape, bbox);

  hmap::export_banner_png("ex_path_smooth.png",
                          {z1, z2, z3},
                          hmap::Cmap::INFERNO);
}
