#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {2.f, 2.f};
  int        seed = 0;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2, shape, kw, seed);
  hmap::remap(z);

  // --- GPU phase field

  float       kp_global = 16.f;
  hmap::Array phi2 = hmap::gpu::phase_field(z, seed, kp_global);
  hmap::remap(phi2);

  hmap::export_banner_png("ex_phase_field.png", {z, phi2}, hmap::Cmap::JET);
}
