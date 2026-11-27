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
  path.resample_uniform();

  auto z1 = path.to_array(shape, bbox);

  hmap::Path path_c;

  //
  float ratio = 0.2f;

  path_c = hmap::meanderize(path, ratio);
  auto z2 = path_c.to_array(shape, bbox);

  //
  ratio = 0.3f;
  float noise_ratio = 0.1f;
  int   iterations = 3;

  path.set_closed(true);
  path_c = hmap::meanderize(path, ratio, noise_ratio, seed, iterations);
  auto z3 = path_c.to_array(shape, bbox);

  hmap::export_banner_png("ex_path_meanderize.png",
                          {z1, z2, z3},
                          hmap::Cmap::INFERNO);
}
