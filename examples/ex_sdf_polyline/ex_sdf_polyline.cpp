#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array noise = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      kw,
                                      seed,
                                      8,
                                      0.f);
  hmap::remap(noise, 0.f, 1.f);

  // base path
  hmap::BBox2D bbox = {1.f, 2.f, -0.5f, 0.5f};
  hmap::Path path = hmap::Path(9, seed, hmap::BBox2D(1.2f, 1.8f, -0.3f, 0.3f));
  path.reorder_nns();
  // path.fractalize(2, seed);
  path.set_values_from_array(noise, bbox);

  auto zp = hmap::Array(shape);
  path.to_array(zp, bbox);

  hmap::UvRect uv_rect = {1.f, 2.f, -0.5f, 0.5f};
  hmap::Array  sdf_2d = hmap::sdf_2d_polyline(path, shape, uv_rect);
  hmap::remap(sdf_2d);

  hmap::Array sdf_2d_bezier = hmap::sdf_2d_polyline_bezier(path,
                                                           shape,
                                                           uv_rect);
  hmap::remap(sdf_2d_bezier);

  hmap::export_banner_png("ex_sdf_polyline.png",
                          {zp, sdf_2d, sdf_2d_bezier},
                          hmap::Cmap::INFERNO);
}
