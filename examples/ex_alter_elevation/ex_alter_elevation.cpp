#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 1;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);

  // modification point (x, y) in [0, 1]x[0, 1] and a value
  // corresponding to the relative elevation modification
  hmap::Cloud cloud = hmap::Cloud();
  cloud.push_back(hmap::Point(0.2f, 0.5f, -1.f));
  cloud.push_back(hmap::Point(0.6f, 0.2f, 1.f));

  hmap::Array z1 = z0;
  hmap::alter_elevation(z1, cloud, 32, 2.f);

  z1.to_png("out.png", hmap::Cmap::INFERNO);

  z0.infos();
  z1.infos();

  hmap::export_banner_png("ex_alter_elevation.png",
                          {z0, z1},
                          hmap::Cmap::INFERNO);
}
