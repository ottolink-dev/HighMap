#include "highmap.hpp"

int main(void)
{
  hmap::init_openmp();

  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};
  glm::vec2 kw = {2.f, 2.f};
  int       seed = 0;

  hmap::Array z1 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z1);

  // --- create path

  hmap::Path path = hmap::find_cut_path_midpoint(
      z1,
      hmap::DomainBoundary::BOUNDARY_LEFT,
      hmap::DomainBoundary::BOUNDARY_RIGHT,
      seed);

  path = hmap::decimate_vw(path, 30);
  path = hmap::bspline(path, 50);

  hmap::Array z2 = hmap::Array(shape);
  path.to_array(z2);

  // path.smooth(10);

  // hmap::Array z2 = hmap::Array(shape);
  path.to_array(z2);

  // --- trenchs

  auto z3 = z1;

  hmap::Array dr = 1.f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                         shape,
                                         {8.f, 8.f},
                                         ++seed,
                                         8,
                                         0.f);

  float width = 0.05f;
  hmap::trench(z3, path, width, true, false);

  z3.dump();

  hmap::export_banner_png("ex_trench.png", {z1, z2, z3}, hmap::Cmap::JET, true);
}
