#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);

  hmap::Tensor tex = hmap::colorize_slope_height_heatmap(z,
                                                         hmap::Cmap::VIRIDIS);
  hmap::Array  r = tex.get_slice(0);
  hmap::Array  g = tex.get_slice(1);
  hmap::Array  b = tex.get_slice(2);

  hmap::ColorAdjust params;
  params.contrast = 1.5f;
  params.saturation = 1.2f;
  hmap::color_adjust(r, g, b, params, {0, 0});

  hmap::Tensor tex_adj(z.shape, 3);
  tex_adj.set_slice(0, r);
  tex_adj.set_slice(1, g);
  tex_adj.set_slice(2, b);
  tex_adj.to_png("ex_color_adjust.png");
}
