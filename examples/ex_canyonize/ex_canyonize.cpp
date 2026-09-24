#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {3.f, 3.f};
  int        seed = 42;

  // Base terrain with ridges
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed++);
  hmap::remap(z, 0.f, 1.f);

  // Optional high-frequency noise for lateral strata distortion
  hmap::Array noise = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      glm::vec2(2.f, 2.f),
                                      seed++);
  hmap::remap(noise, 0.f, 0.1f);

  auto z_canyon = z;
  hmap::canyonize(z_canyon,
                  seed++,
                  12,    // 6 strata levels
                  2.f,   // convex strata height ratio relative to concave
                  2.f,   // convex exponent for base and alternate levels
                  0.5f,  // concave exponent for alternate levels
                  0.15f, // clamp min to flatten valley / river floor
                  0.08f, // smooth clamp transition
                  0.7f,  // strata height fluctuation
                  &noise);

  z_canyon.dump();

  hmap::export_banner_png("ex_canyonize.png",
                          {z, z_canyon},
                          hmap::Cmap::TERRAIN,
                          true);
}
