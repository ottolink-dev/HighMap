#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 42;

  // base terrain using fractional Brownian motion
  hmap::Array z1 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z1);

  // accentuate ridges using gradient-directed resampling
  hmap::Array z2 = z1;
  hmap::Array z3 = z1;

  float strength = 0.01f;
  int   ir = 32;

  hmap::ridge_accentuate(z2, strength, ir);
  hmap::gpu::ridge_accentuate(z3, strength, ir);

  z1.dump("z1.png");
  z2.dump("z2.png");
  z3.dump("z3.png");

  // export before and after banner
  hmap::export_banner_png("ex_ridge_accentuate.png",
                          {z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
