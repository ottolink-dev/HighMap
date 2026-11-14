#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};
  glm::vec2     kw = {2.f, 2.f};
  std::uint32_t seed = 1;

  hmap::Array z = hmap::Array(shape);
  glm::vec4   bbox = hmap::unit_square_bbox();

  float mu = 1.02f;

  hmap::ParberryFunction p(kw, seed, mu);

  hmap::FbmFunction f(p, 8, 0.7f, 0.5f, 2.f);

  fill_array_using_xy_function(z, bbox, nullptr, nullptr, nullptr, nullptr, f);

  z.to_png("out.png", hmap::Cmap::TERRAIN, true);
  z.infos();
}
