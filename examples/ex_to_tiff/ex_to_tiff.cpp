#include "highmap.hpp"

int main(void)
{
  glm::ivec2  shape = {256, 256};
  glm::vec2   kw = {4.f, 4.f};
  hmap::Array z = hmap::noise(hmap::NoiseType::SIMPLEX2, shape, kw, 1);
  z.to_tiff("ex_to_tiff_output.tiff");
  return 0;
}
