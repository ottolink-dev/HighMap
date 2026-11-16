#include <iostream>

#include "highmap.hpp"
#include "highmap/internal/timer.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2    shape = {256, 256};
  glm::vec2     kw = {2.f, 2.f};
  std::uint32_t seed = 0;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);

  hmap::Array mask_fixed_values = z0;
  hmap::clamp_min(mask_fixed_values, 0.f);

  mask_fixed_values.infos();

  // 1. Isotropic harmonic interpolation (Dx = Dy = 1.0)
  hmap::Timer::Start("GPU Isotropic");
  hmap::Array z_iso = hmap::gpu::harmonic_interpolation(z0, mask_fixed_values);
  hmap::Timer::Stop("GPU Isotropic");

  // 2. Anisotropic harmonic interpolation with 2D varying diffusion
  // coefficients (Dx, Dy)
  hmap::Array dx(shape, 100.0f);
  hmap::Array dy(shape, 1.0f);

  hmap::Timer::Start("GPU Anisotropic");
  hmap::Array z_aniso = hmap::gpu::harmonic_interpolation(z0,
                                                          mask_fixed_values,
                                                          dx,
                                                          dy);
  hmap::Timer::Stop("GPU Anisotropic");

  hmap::export_banner_png("ex_harmonic_interpolation.png",
                          {z0, mask_fixed_values, z_iso, z_aniso},
                          hmap::Cmap::JET,
                          false);

  return 0;
}
