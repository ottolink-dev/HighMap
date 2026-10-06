#include "highmap.hpp"

int main(void)
{
  hmap::log::warn(
      "hydraulic_diffusion is deprecated and will be removed at some point.");

  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  auto        z0 = z;

  float c_diffusion = 0.1f;
  int   iterations = 100;

  auto  dn = gradient_norm(z);
  float talus_ref = 1.1f * dn.max();

  hmap::hydraulic_diffusion(z, c_diffusion, talus_ref, iterations);

  hmap::export_banner_png("ex_hydraulic_diffusion.png",
                          {z0, z},
                          hmap::Cmap::TERRAIN,
                          true);
}
