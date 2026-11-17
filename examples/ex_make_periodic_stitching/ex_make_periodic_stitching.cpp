#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  float       overlap = 0.5f;

  auto zp = hmap::make_periodic_stitching(z, overlap);

  // do some tiling to check periodicity
  auto zt = zp;
  zt = hstack(zt, zt);
  zt = vstack(zt, zt);

  zt.infos();

  hmap::export_banner_png("ex_make_periodic_stitching0.png",
                          {z, zp},
                          hmap::Cmap::INFERNO);

  zt.to_png("ex_make_periodic_stitching1.png", hmap::Cmap::INFERNO);

  auto zt_grad = hmap::gradient_norm(zt);
  zt_grad.to_png("ex_make_periodic_stitching_grad.png", hmap::Cmap::INFERNO);
}
