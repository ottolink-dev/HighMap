#include <iostream>

#include "highmap.hpp"

int main(void)
{
  float a = 0.5f;
  float b = 0.7f;
  float val = hmap::maximum_smooth(a, b, 0.2f);

  // Also generate an image for Doxygen
  glm::ivec2  shape = {256, 256};
  glm::vec2   kw = {4.f, 4.f};
  hmap::Array z = hmap::noise(hmap::NoiseType::SIMPLEX2, shape, kw, 1);
  hmap::remap(z);

  hmap::Array z_smooth(shape);
  for (int i = 0; i < shape.x; ++i)
  {
    for (int j = 0; j < shape.y; ++j)
    {
      z_smooth(i, j) = hmap::maximum_smooth(z(i, j), 0.5f, 0.2f);
    }
  }

  hmap::export_banner_png("ex_maximum_smooth_scalar.png",
                          {z, z_smooth},
                          hmap::Cmap::VIRIDIS);
  return 0;
}
