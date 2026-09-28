#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 42;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw,
                                   seed,
                                   8,
                                   0.7f);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 1.f);
  hmap::remap(z0);

  int nparticles = int(1.f * shape.x * shape.y);
  int iterations = 1;

  // 1. Without trail (standard hydraulic_particle)
  auto z1 = z0;
  hmap::gpu::hydraulic_particle(z1,
                                nparticles,
                                seed,
                                /* p_bedrock */ nullptr,
                                /* p_moisture_map */ nullptr,
                                /* p_elevation_shift */ nullptr,
                                /* p_erosion_map */ nullptr,
                                /* p_deposition_map */ nullptr,
                                /* c_capacity */ 10.f,
                                /* c_erosion */ 0.05f,
                                /* c_deposition */ 0.05f,
                                /* c_inertia */ 0.01f,
                                /* c_gravity */ 1.f,
                                /* drag_rate */ 0.001f,
                                /* evap_rate */ 0.001f,
                                /* talus_slope */ 2.f,
                                /* collapse_rate */ 0.1f,
                                /* iterations */ iterations);

  // 2. With trail (moderate trail attraction)
  auto        z2 = z0;
  hmap::Array trail2(shape, 0.f);
  hmap::gpu::hydraulic_particle_trail(z2,
                                      nparticles,
                                      seed,
                                      /* p_bedrock */ nullptr,
                                      /* p_moisture_map */ nullptr,
                                      /* p_elevation_shift */ nullptr,
                                      /* p_erosion_map */ nullptr,
                                      /* p_deposition_map */ nullptr,
                                      /* p_trail_map */ &trail2,
                                      /* c_capacity */ 10.f,
                                      /* c_erosion */ 0.05f,
                                      /* c_deposition */ 0.05f,
                                      /* c_inertia */ 0.01f,
                                      /* c_gravity */ 1.f,
                                      /* drag_rate */ 0.001f,
                                      /* evap_rate */ 0.001f,
                                      /* talus_slope */ 2.f,
                                      /* collapse_rate */ 0.1f,
                                      /* c_trail_deposit */ 1.f,
                                      /* c_trail_attraction */ 0.4f,
                                      /* trail_evap_rate */ 0.1f,
                                      /* iterations */ iterations);

  // 3. With trail (strong trail attraction)
  auto        z3 = z0;
  hmap::Array trail3(shape, 0.f);
  hmap::gpu::hydraulic_particle_trail(z3,
                                      nparticles,
                                      seed,
                                      /* p_bedrock */ nullptr,
                                      /* p_moisture_map */ nullptr,
                                      /* p_elevation_shift */ nullptr,
                                      /* p_erosion_map */ nullptr,
                                      /* p_deposition_map */ nullptr,
                                      /* p_trail_map */ &trail3,
                                      /* c_capacity */ 10.f,
                                      /* c_erosion */ 0.05f,
                                      /* c_deposition */ 0.05f,
                                      /* c_inertia */ 0.01f,
                                      /* c_gravity */ 1.f,
                                      /* drag_rate */ 0.001f,
                                      /* evap_rate */ 0.001f,
                                      /* talus_slope */ 2.f,
                                      /* collapse_rate */ 0.1f,
                                      /* c_trail_deposit */ 1.f,
                                      /* c_trail_attraction */ 0.7f,
                                      /* trail_evap_rate */ 0.5f,
                                      /* iterations */ iterations);

  // Remap trail for visualization
  trail3.infos();

  hmap::Array trail_vis = trail3;
  hmap::remap(trail_vis);

  // --- Output

  z0.dump("out0.png");
  z1.dump("out1.png");
  z2.dump("out2.png");
  z3.dump("out3.png");
  trail_vis.dump("out_trail.png");

  hmap::export_banner_png("ex_hydraulic_particle_trail.png",
                          {z0, z1, z2, z3, trail_vis},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
