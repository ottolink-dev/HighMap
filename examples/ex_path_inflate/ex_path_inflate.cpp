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

  path = hmap::fractalize(path, 2, seed, 0.2f);
  path.resample_interp(200);

  auto z1 = path.to_array(shape, bbox);

  // --- open

  hmap::Array z2(shape);

  for (const auto &radius : {0.f, 0.1f, 0.2f, 0.3f})
  {
    hmap::Path p = hmap::inflate(path, radius);
    p.set_values(radius + 0.7f);
    p.to_array(z2, bbox);
  }

  // --- closed

  path.set_closed(true);

  hmap::Array z3(shape);

  for (const auto &radius : {0.f, 0.1f, 0.2f, 0.3f})
  {
    hmap::Path p = hmap::inflate(path, radius);
    p.set_values(radius + 0.7f);
    p.to_array(z3, bbox);
  }

  hmap::export_banner_png("ex_path_inflate.png",
                          {z1, z2, z3},
                          hmap::Cmap::INFERNO);
}
