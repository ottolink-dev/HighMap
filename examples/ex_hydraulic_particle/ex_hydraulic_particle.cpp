#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};
  glm::vec2 kw = {4.f, 4.f};
  int       seed = 0;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw,
                                   seed,
                                   8,
                                   0.f);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 1.f);
  hmap::remap(z0);

  int nparticles = int(0.5f * shape.x * shape.y);

  auto z1 = z0;
  hmap::gpu::hydraulic_particle(z1,
                                nparticles,
                                seed,
                                /* p_bedrock */ nullptr,
                                /* p_erodibility */ nullptr,
                                /* p_moisture_map */ nullptr,
                                /* p_elevation_shift */ nullptr,
                                /* p_erosion_map */ nullptr,
                                /* p_deposition_map */ nullptr,
                                /* c_capacity */ 10.f,
                                /* c_erosion */ 0.1f,
                                /* c_deposition */ 0.1f,
                                /* c_gravity */ 1.f,
                                /* drag_rate */ 0.001f,
                                /* evap_rate */ 0.001f,
                                /* talus_slope */ 2.f,
                                /* collapse_rate */ 0.25f,
                                /* iterations */ 1);

  z1.dump();

  auto z2 = z0;
  auto zb2 = hmap::generate_bedrock(z2,
                                    /* elevation_strength */ 0.05f,
                                    /* slope_strength */ 0.f,
                                    /* slope_limit */ 0.f);
  hmap::gpu::hydraulic_particle(z2, nparticles, seed, &zb2);

  auto z3 = z0;
  auto zb3 = hmap::generate_bedrock(z3,
                                    /* elevation_strength */ 0.f,
                                    /* slope_strength */ 0.1f,
                                    /* slope_limit */ 2.f / shape.x);
  hmap::gpu::hydraulic_particle(z3, nparticles, seed, &zb3);

  auto z4 = z0;
  hmap::gpu::hydraulic_particle_multiscale(
      z4,
      seed,
      /* steps_per_level */ {16, 8, 4, 2, 1},
      /* p_bedrock */ nullptr,
      /* p_erodibility */ nullptr,
      /* p_moisture_map */ nullptr,
      /* p_elevation_shift */ nullptr,
      /* p_erosion_map */ nullptr,
      /* p_deposition_map */ nullptr,
      /* particles_ratio */ 0.5f,
      /* c_capacity */ 10.f,
      /* c_erosion */ 0.05f,
      /* c_deposition */ 0.05f,
      /* c_inertia */ 0.0f,
      /* c_gravity */ 1.f,
      /* drag_rate */ 0.001f,
      /* evap_rate */ 0.001f,
      /* talus_slope */ 2.f,
      /* collapse_rate */ 0.25f,
      /* mix */ 0.8f);

  // --- output

  z0.dump("out0.png");
  z1.dump("out1.png");
  z2.dump("out2.png");
  z3.dump("out3.png");
  z4.dump("out4.png");

  hmap::export_banner_png("ex_hydraulic_particle.png",
                          {z0, z1, z2, z3, z4},
                          hmap::Cmap::TERRAIN,
                          true);
}
