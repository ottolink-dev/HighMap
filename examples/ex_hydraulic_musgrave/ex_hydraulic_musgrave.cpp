#include "highmap.hpp"

int main(void)
{
  hmap::log::warn(
      "hydraulic_musgrave is deprecated and will be removed at some point.");

  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  auto        z0 = z;

  hmap::hydraulic_musgrave(z);

  hmap::export_banner_png("ex_hydraulic_musgrave.png",
                          {z0, z},
                          hmap::Cmap::TERRAIN,
                          true);
}
