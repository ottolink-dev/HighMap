#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 6;
  int        npoints = 10;

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path   path = hmap::Path(npoints,
                               seed,
                               hmap::BBox2D(0.05f, 0.95f, 0.05f, 0.95f));
  path.reorder_nns();
  path.set_closed(true);

  path = hmap::fractalize(path, 2, seed, 0.2f);
  path.resample_interp(200);

  auto z1 = path.to_array(shape, bbox);

  // --- Concentric uniform scale (shrinking creates margin from boundary)

  hmap::Array z2(shape);
  for (float s : {1.0f, 0.75f, 0.5f, 0.25f})
  {
    hmap::Path p = hmap::scale(path, s);
    p.to_array(z2, bbox);
  }

  // --- Non-uniform scale (independent X and Y scaling)

  hmap::Array z3(shape);
  for (float s : {1.0f, 0.75f, 0.5f, 0.25f})
  {
    hmap::Path p = hmap::scale(path, glm::vec2(s, s * 0.5f));
    p.to_array(z3, bbox);
  }

  hmap::export_banner_png("ex_path_scale.png",
                          {z1, z2, z3},
                          hmap::Cmap::INFERNO);

  return 0;
}
