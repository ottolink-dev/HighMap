#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 2;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  kw,
                                  seed,
                                  8,
                                  0.7f);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 1.f);
  hmap::remap(z);

  auto map =
      hmap::disk(shape, 0.05f, 32.f, nullptr, nullptr, nullptr, {0.5f, 0.45f});

  float amount = 1e-1f;
  int   iterations = 500000;

  hmap::Array water_depth = hmap::gpu::flow_simulation_viscous(z,
                                                               amount,
                                                               map,
                                                               iterations);

  z.dump("z.png");
  water_depth.dump("depth.png");

  hmap::export_banner_png("ex_flow_simulation_viscous.png",
                          {z, z + water_depth, water_depth},
                          hmap::Cmap::JET);
}
