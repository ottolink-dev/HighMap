#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 3;
  int        npoints = 30;

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path   path = hmap::Path(npoints,
                               seed,
                               hmap::BBox2D(0.1f, 0.9f, 0.1f, 0.9f));
  path.reorder_nns();
  path = hmap::fractalize(path, 1, seed);

  int ntarget = 15;

  auto z1 = hmap::Array(shape);
  path.to_array(z1, bbox);

  // Visvalingam-Whyatt algo
  auto path2 = hmap::decimate_vw(path, ntarget);
  auto z2 = path2.to_array(shape, bbox);

  hmap::assert_start_end_points(path, path2);

  path.set_closed(true);
  auto path3 = hmap::decimate_vw(path, ntarget);
  auto z3 = path3.to_array(shape, bbox);

  hmap::export_banner_png("ex_path_decimate.png",
                          {z1, z2, z3},
                          hmap::Cmap::INFERNO);
}
