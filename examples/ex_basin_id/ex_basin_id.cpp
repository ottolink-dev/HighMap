#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {8.f, 8.f};
  int        seed = 0;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CONE, 1.f);

  hmap::Array id1 = hmap::basin_id(z, false);
  hmap::remap(id1);

  hmap::Array id2 = hmap::basin_id(z, true);
  hmap::remap(id2);

  hmap::export_banner_png("ex_basin_id.png", {z, id1, id2}, hmap::Cmap::JET);
}
