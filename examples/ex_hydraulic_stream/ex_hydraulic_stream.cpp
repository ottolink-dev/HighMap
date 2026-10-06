#include "highmap.hpp"

int main(void)
{
  hmap::log::warn(
      "hydraulic_stream is deprecated and will be removed at some point.");

  hmap::gpu::init_opencl();
  hmap::init_openmp();

  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};
  glm::vec2 kw = {2.f, 2.f};
  int       seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);
  auto z0 = z;

  float       c_erosion = 0.1f;
  float       talus_ref = 5.f / (float)shape.x;
  int         iradius = 64;
  hmap::Array z_bedrock = hmap::local_min(z, iradius);

  auto z1 = z;
  hmap::hydraulic_stream(z1, c_erosion, talus_ref);

  auto z2 = z;
  int  ir = 5;

  auto erosion_map = hmap::Array(shape);
  auto moisture_map = z;

  hmap::hydraulic_stream(z2,
                         c_erosion,
                         talus_ref,
                         &z_bedrock,
                         &moisture_map,
                         &erosion_map,
                         ir);

  // log scale
  auto z3 = z;
  int  deposition_ir = 32;
  c_erosion = 0.2f;
  hmap::gpu::hydraulic_stream_log(z3, c_erosion, talus_ref, deposition_ir);

  z3.dump();

  // multiscale log scale
  auto z4 = z;
  hmap::gpu::hydraulic_stream_log_multiscale(z4,
                                             c_erosion,
                                             talus_ref,
                                             /* levels */ 3,
                                             deposition_ir);

  z4.dump();

  hmap::export_banner_png("ex_hydraulic_stream0.png",
                          {z0, z1, z2, z3, z4},
                          hmap::Cmap::TERRAIN,
                          true);

  erosion_map.to_png("ex_hydraulic_stream1.png", hmap::Cmap::INFERNO);
}
