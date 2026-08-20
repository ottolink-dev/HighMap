#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);

  hmap::Array z1 = z;
  hmap::chop(z1, 0.4f);

  hmap::export_banner_png("ex_chop.png", {z, z1}, hmap::Cmap::VIRIDIS);
}
