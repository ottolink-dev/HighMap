#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};

  hmap::Array z1 = hmap::noise(hmap::NoiseType::SIMPLEX2, shape, kw, 1);
  hmap::Array z2 = hmap::noise(hmap::NoiseType::SIMPLEX2, shape, kw, 2);
  hmap::remap(z1);
  hmap::remap(z2);

  hmap::Array z_smooth = hmap::maximum_smooth(z1, z2, 0.2f);

  hmap::export_banner_png("ex_maximum_smooth.png",
                          {z1, z2, z_smooth},
                          hmap::Cmap::VIRIDIS);
}
