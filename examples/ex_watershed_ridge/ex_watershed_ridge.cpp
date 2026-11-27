#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 0;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed++);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 2.f);
  hmap::remap(z0);

  float       noise_amp = 0.03f;
  hmap::Array dx = noise_amp * hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                               shape,
                                               4.f * kw,
                                               seed);

  // scale width and amplitude (smaller at the domain center)
  hmap::Array scaling = 1.f - hmap::cubic_pulse(shape);

  float amplitude = 0.2f;
  float width = 24.f;
  float exponent = 0.5f;
  int   prefilter_ir = 4;

  auto z1 = hmap::gpu::watershed_ridge(z0,
                                       amplitude,
                                       width,
                                       exponent,
                                       prefilter_ir,
                                       false,
                                       &dx,
                                       &dx,
                                       &scaling);

  auto z2 = hmap::gpu::watershed_ridge(z0,
                                       amplitude,
                                       width,
                                       exponent,
                                       prefilter_ir,
                                       true,
                                       &dx,
                                       &dx,
                                       &scaling);

  hmap::export_banner_png("ex_watershed_ridge.png",
                          {z0, z1, z2},
                          hmap::Cmap::TERRAIN,
                          true);
}
