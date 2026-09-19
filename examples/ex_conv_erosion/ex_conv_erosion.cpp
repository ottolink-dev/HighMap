#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 4;

  shape = {512, 512};
  // shape = {1024, 1024};

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 0.5f);

  auto z1 = z0;

  int   iterations = 20;
  int   particle_count = 2500;
  int   ir_min = 1;
  int   ir_max = 1;
  float size_distrib_exp = 2.f;
  float erosion_strength = 0.01f;
  float randomness = 0.002f;
  float exit_forcing = 0.01f;

  hmap::gpu::conv_erosion(z1,
                          seed,
                          iterations,
                          particle_count,
                          ir_min,
                          ir_max,
                          size_distrib_exp,
                          erosion_strength,
                          randomness,
                          exit_forcing);

  // hmap::gpu::deposition_fill_holes(z1, 4, 1.f);

  z1.dump("out.png");

  hmap::export_banner_png("ex_conv_erosion.png",
                          {z0, z1},
                          hmap::Cmap::TERRAIN,
                          true);
}
