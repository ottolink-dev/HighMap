#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;
  glm::vec4  bbox = hmap::unit_square_bbox();

  hmap::PerlinFunction p = hmap::PerlinFunction(kw, seed);

  auto n = 0.4f * hmap::noise(hmap::NoiseType::PERLIN, shape, kw, ++seed);

  auto z1 = hmap::Array(shape);
  auto z2 = hmap::Array(shape);

  fill_array_using_xy_function(z1, bbox, nullptr, &n, &n, nullptr, p);

  int subsampling = 8;

  fill_array_using_xy_function(z2,
                               bbox,
                               nullptr, // ctrl_param
                               &n,
                               &n,
                               nullptr, // stretching
                               p,
                               subsampling);

  hmap::export_banner_png("ex_fill_array_using_xy_function.png",
                          {z1, z2},
                          hmap::Cmap::JET);
}
