#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 42;

  // Generate base terrain with noise and erosion-like smoothing
  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z0, 0.f, 1.f);

  // Apply adaptive relief enhancement filter on CPU
  hmap::Array z_cpu = z0;
  hmap::adaptive_relief(z_cpu, 1.5f, 0.0f);

  // Apply adaptive relief enhancement filter on GPU (OpenCL)
  hmap::Array z_gpu = z0;
  hmap::gpu::adaptive_relief(z_gpu, 1.5f, 0.0f);

  // Stack original, CPU enhanced, and GPU enhanced side-by-side
  hmap::export_banner_png("ex_adaptive_relief.png",
                          {z0, z_cpu, z_gpu},
                          hmap::Cmap::JET);

  return 0;
}
