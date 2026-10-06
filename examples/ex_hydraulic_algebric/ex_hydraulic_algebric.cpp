#include "highmap.hpp"

int main(void)
{
  hmap::log::warn(
      "hydraulic_algebric is deprecated and will be removed at some point.");

  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  auto        z0 = z;

  auto  talus = gradient_talus(z);
  float talus_ref = 0.2f * talus.max();
  float c_erosion = 0.05f;
  float c_deposition = 0.05f;
  int   iterations = 1;
  int   ir = 16;

  hmap::hydraulic_algebric(z,
                           talus_ref,
                           ir,
                           nullptr,
                           nullptr,
                           nullptr,
                           c_erosion,
                           c_deposition,
                           iterations);

  hmap::export_banner_png("ex_hydraulic_algebric.png",
                          {z0, z},
                          hmap::Cmap::TERRAIN,
                          true);
}
