#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  const glm::ivec2 shape = {512, 512};
  const glm::vec2  kw_base = {2.f, 2.f};
  const int        seed = 42;

  // Base terrain heightmap
  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw_base,
                                   seed);
  hmap::remap(z0, 0.f, 1.f);

  const float     kw = 8.f;
  const float     kw_fine = 16.f;
  const glm::vec2 jitter = {0.6f, 0.6f};
  const float     factor_mild = 0.8f;
  const float     factor_strong = 1.5f;
  const float     shape_factor_none = 0.f;
  const float     shape_factor_smooth = 0.5f;
  const float     shape_factor_accentuated = 1.5f;
  const float     mask_sigma = 0.25f;

  // Jagged filter with hard cell boundaries (shape_factor = 0.0)
  auto z1 = hmap::gpu::jagged(z0,
                              kw,
                              seed,
                              jitter,
                              factor_mild,
                              shape_factor_none);

  // Jagged filter with edge distance falloff (shape_factor = 1.0)
  auto z2 = hmap::gpu::jagged(z0,
                              kw,
                              seed,
                              jitter,
                              factor_mild,
                              shape_factor_smooth);

  z2.dump();

  // Jagged filter with higher frequency and mask
  hmap::Array mask = hmap::gaussian_pulse(shape, mask_sigma);
  auto        z3 = hmap::gpu::jagged(z0,
                              kw_fine,
                              seed,
                              jitter,
                              factor_strong,
                              shape_factor_accentuated,
                              &mask);

  hmap::export_banner_png("ex_jagged.png",
                          {z0, z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
