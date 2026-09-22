#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {1024, 1024};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 0;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);

  float talus = 3.f / shape.x;
  int   ir = 0;
  float amplitude = 2.f;
  float gain = 2.f;

  // isotropic cliff exaggeration
  auto z_cliff = z;
  hmap::recast_cliff(z_cliff, talus, ir, amplitude, gain);
  hmap::remap(z_cliff);

  // directional cliff exaggeration (angle = 45 deg)
  auto  z_cliff_dir = z;
  float angle = 45.f;

  hmap::recast_cliff_directional(z_cliff_dir,
                                 talus,
                                 ir,
                                 amplitude,
                                 angle,
                                 gain);
  hmap::remap(z_cliff_dir);

  // directional cliff exaggeration with variable (local) angle
  auto        z_cliff_var = z;
  hmap::Array angle_field = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                            shape,
                                            glm::vec2(2.f, 2.f),
                                            seed + 1);
  hmap::remap(angle_field, 0.f, 360.f);

  hmap::recast_cliff_directional(z_cliff_var,
                                 talus,
                                 ir,
                                 amplitude,
                                 angle_field,
                                 gain);
  hmap::remap(z_cliff_var);

  hmap::export_banner_png("ex_recast_cliff.png",
                          {z, z_cliff, z_cliff_dir, z_cliff_var},
                          hmap::Cmap::TERRAIN,
                          true);
}
