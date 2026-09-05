#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl(); // for hydraulic erosion

  glm::ivec2 shape_hr = {512, 512};
  int        factor = 4;
  int        seed = 1;

  // high-resolution exemplar providing the fine details
  hmap::Array exemplar = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                         shape_hr,
                                         {2.f, 2.f},
                                         seed + 1);

  hmap::gpu::hydraulic_particle(exemplar, 500000, seed);
  hmap::remap(exemplar);

  // smooth, low-resolution input terrain to amplify
  glm::ivec2  shape_lr = {shape_hr.x / factor, shape_hr.y / factor};
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape_lr,
                                  {4.f, 4.f},
                                  ++seed,
                                  2);
  hmap::remap(z);

  int patch_size = 16;
  int analysis_stride = patch_size / 8;
  int synthesis_stride = patch_size / 2;

  hmap::Array zsr = hmap::terrain_super_resolution(z,
                                                   exemplar,
                                                   factor,
                                                   patch_size,
                                                   analysis_stride,
                                                   synthesis_stride);

  zsr.dump();

  // plain upsampling of the input, for comparison
  hmap::Array zup = z.resample_to_shape_bicubic(shape_hr);

  hmap::export_banner_png("ex_terrain_super_resolution.png",
                          {exemplar, zup, zsr},
                          hmap::Cmap::TERRAIN,
                          true);
}
