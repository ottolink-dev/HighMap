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

  const float     amp = 0.15f;
  const int       octaves = 4;
  const float     persistence = 0.5f;
  const float     lacunarity = 2.0f;
  const glm::vec2 jitter = {0.6f, 0.6f};

  // 1. Isotropic jagged fBm (kw = 6)
  auto z1 = hmap::gpu::jagged_fbm(z0,
                                  6.f,
                                  amp,
                                  seed,
                                  octaves,
                                  persistence,
                                  lacunarity,
                                  /* switch_kx_ky */ false,
                                  jitter);

  // 2. Anisotropic jagged fBm without switching (kw = {12, 3})
  const glm::vec2 kw_aniso = {12.f, 3.f};
  auto            z2 = hmap::gpu::jagged_fbm(z0,
                                  kw_aniso,
                                  amp,
                                  seed,
                                  octaves,
                                  persistence,
                                  lacunarity,
                                  /* switch_kx_ky */ false,
                                  jitter);

  // 3. Anisotropic jagged fBm with switching kx and ky at each octave
  auto z3 = hmap::gpu::jagged_fbm(z0,
                                  kw_aniso,
                                  amp,
                                  seed,
                                  octaves,
                                  persistence,
                                  lacunarity,
                                  /* switch_kx_ky */ true,
                                  jitter);

  hmap::export_banner_png("ex_jagged_fbm.png",
                          {z0, z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
