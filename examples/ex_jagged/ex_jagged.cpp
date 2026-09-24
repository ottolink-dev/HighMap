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
  const float     amp_mild = 0.12f;
  const float     amp_strong = 0.22f;
  const float     gamma_plateau = 0.5f;
  const float     gamma_cusp = 2.0f;
  const float     mask_sigma = 0.25f;

  // Jagged bubble domes with plateau profile (gamma = 0.5)
  auto z1 = hmap::gpu::jagged(z0, kw, amp_mild, seed, jitter, gamma_plateau);

  // Jagged bubble domes with standard parabolic profile (gamma = 1.0)
  auto z2 = hmap::gpu::jagged(z0, kw, amp_mild, seed, jitter, 1.0f);

  // Jagged bubble domes with cusp profile and mask (gamma = 2.0)
  hmap::Array mask = hmap::gaussian_pulse(shape, mask_sigma);
  auto        z3 = hmap::gpu::jagged(z0,
                              kw_fine,
                              amp_strong,
                              seed,
                              jitter,
                              gamma_cusp,
                              30.f,
                              &mask);

  hmap::export_banner_png("ex_jagged.png",
                          {z0, z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
