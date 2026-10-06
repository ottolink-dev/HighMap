#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {2.f, 2.f};
  int        seed = 2;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z0);
  auto z1 = z0;

  int         iterations = 500;
  hmap::Array talus(shape, 4.f / (float)shape.x); // for thermal

  hmap::gpu::hydraulic_schott(z1,
                              iterations,
                              talus,
                              0.5f,
                              0.1f,
                              0.05f,
                              0.8f,
                              0.8f,
                              1.3f,
                              1.f,
                              1.f);

  // erosion only
  auto z2 = z0;
  hmap::gpu::hydraulic_schott_erosion(z2, 20);

  hmap::export_banner_png("ex_hydraulic_schott.png",
                          {z0, z1, z2},
                          hmap::Cmap::TERRAIN,
                          true);
}
