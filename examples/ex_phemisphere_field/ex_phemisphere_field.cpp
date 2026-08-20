#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z = hmap::gpu::hemisphere_field(shape, kw, seed);
  hmap::remap(z);

  z.to_png("ex_phemisphere_field.png", hmap::Cmap::VIRIDIS);
}
