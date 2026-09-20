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
  const float     gamma_mild = 0.5f;
  const float     gamma_strong = 1.5f;
  const float     shape_gamma_none = 0.f;
  const float     shape_gamma_smooth = 0.5f;
  const float     shape_gamma_accentuated = 1.5f;
  const float     mask_sigma = 0.25f;

  // Jagged filter with hard cell boundaries (shape_gamma = 0.0)
  auto z1 = hmap::gpu::jagged(z0,
                              kw,
                              seed,
                              jitter,
                              gamma_mild,
                              shape_gamma_none);

  // Jagged filter with edge distance falloff (shape_gamma = 0.5)
  auto z2 = hmap::gpu::jagged(z0,
                              kw,
                              seed,
                              jitter,
                              gamma_mild,
                              shape_gamma_smooth);

  // Jagged filter with higher frequency, angle, and mask
  hmap::Array mask = hmap::gaussian_pulse(shape, mask_sigma);
  auto        z3 = hmap::gpu::jagged(z0,
                              kw_fine,
                              seed,
                              jitter,
                              gamma_strong,
                              shape_gamma_accentuated,
                              1.f,
                              30.f,
                              &mask);

  hmap::export_banner_png("ex_jagged.png",
                          {z0, z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
