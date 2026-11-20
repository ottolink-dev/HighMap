#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {1024, 1024};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 0;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);
  hmap::hydraulic_mise(z);

  int        cubemap_resolution = shape.x / 4;
  float      overlap = 0.25f;
  int        ir = 16;
  hmap::Cmap cmap = hmap::Cmap::TERRAIN;
  bool       splitted = false;

  hmap::export_as_cubemap("ex_export_as_cubemap.png",
                          z,
                          cubemap_resolution,
                          overlap,
                          ir,
                          cmap,
                          splitted);
}
