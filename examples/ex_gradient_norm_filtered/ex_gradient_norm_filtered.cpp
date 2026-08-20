#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);

  hmap::Array g = hmap::gradient_norm_filtered(z, 5);
  hmap::remap(g);

  hmap::export_banner_png("ex_gradient_norm_filtered.png",
                          {z, g},
                          hmap::Cmap::VIRIDIS);
}
