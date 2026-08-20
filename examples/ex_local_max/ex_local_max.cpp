#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);

  hmap::Array z_max = hmap::local_max(z, 5);
  hmap::Array z_min = hmap::local_min(z, 5);

  hmap::export_banner_png("ex_local_max.png",
                          {z, z_max, z_min},
                          hmap::Cmap::VIRIDIS);
}
