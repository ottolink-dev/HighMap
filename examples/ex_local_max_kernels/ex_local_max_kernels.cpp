#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 42;

  // input heightmap
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);

  const int radius = 12;

  // compare exact disk kernel vs separable octagonal approximation vs separable
  // square kernel
  hmap::Array z_disk = hmap::gpu::local_max(z, radius);
  hmap::Array z_octagon = hmap::gpu::local_max_octagon(z, radius);
  hmap::Array z_square = hmap::gpu::local_max_square(z, radius);

  hmap::export_banner_png("ex_local_max_kernels.png",
                          {z, z_disk, z_octagon, z_square},
                          hmap::Cmap::VIRIDIS);
}
