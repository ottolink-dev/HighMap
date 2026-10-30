#include <iostream>

#include "highmap.hpp"

int main(void)
{
  // initialize OpenCL and OpenMP
  hmap::gpu::init_opencl();
  hmap::init_openmp();

  const glm::ivec2 shape = {512, 512};
  const glm::vec2  kw = {3.f, 3.f};
  const int        seed = 42;

  // --- Initial heightfield generation

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw,
                                   seed,
                                   8,
                                   0.f);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 0.8f);
  hmap::remap(z0);

  // --- 1. Hydraulic Particle Erosion (Lagrangian Droplets)

  auto z_particle = z0;
  const int nparticles = static_cast<int>(0.6f * shape.x * shape.y);

  hmap::gpu::hydraulic_particle(
      z_particle,
      nparticles,
      seed,
      /* p_bedrock */ nullptr,
      /* p_moisture_map */ nullptr,
      /* p_elevation_shift */ nullptr,
      /* p_erosion_map */ nullptr,
      /* p_deposition_map */ nullptr,
      /* c_capacity */ 10.f,
      /* c_erosion */ 0.08f,
      /* c_deposition */ 0.08f,
      /* c_inertia */ 0.05f,
      /* c_gravity */ 1.f,
      /* drag_rate */ 0.001f,
      /* evap_rate */ 0.001f,
      /* talus_slope */ 2.f,
      /* collapse_rate */ 0.2f,
      /* iterations */ 4);

  // --- 2. Hydraulic Stream Log Erosion (Eulerian D-inf Stream Power)

  auto z_stream = z0;
  const float c_erosion_stream = 0.2f;
  const float talus_ref = 4.f / static_cast<float>(shape.x);
  const int   deposition_ir = 32;

  hmap::gpu::hydraulic_stream_log(z_stream,
                                  c_erosion_stream,
                                  talus_ref,
                                  deposition_ir,
                                  /* deposition_scale_ratio */ 0.8f);

  // --- 3. Hydraulic MISE (Implicit Stream Power & Sediment Conservation)

  auto z_mise = z0;
  hmap::MiseParams mise_params;
  mise_params.strength = 0.6f;
  mise_params.area_exp = 0.45f;
  mise_params.talus = 3.2f;
  mise_params.debris = 150.f;
  mise_params.deposition_rate = 1.f;
  mise_params.base_res = 128;
  mise_params.seed = seed;

  hmap::hydraulic_mise(z_mise, mise_params);

  // --- Compute erosion difference maps for comparison

  hmap::Array diff_particle = z0 - z_particle;
  hmap::clamp_min(diff_particle, 0.f);

  hmap::Array diff_stream = z0 - z_stream;
  hmap::clamp_min(diff_stream, 0.f);

  hmap::Array diff_mise = z0 - z_mise;
  hmap::clamp_min(diff_mise, 0.f);

  // --- Export individual images and banners with hillshading enabled

  z0.to_png("docs/images/hydraulic_comparison_original.png",
            hmap::Cmap::TERRAIN,
            /* hillshading */ true);
  z_particle.to_png("docs/images/hydraulic_comparison_particle.png",
                    hmap::Cmap::TERRAIN,
                    /* hillshading */ true);
  z_stream.to_png("docs/images/hydraulic_comparison_stream_log.png",
                  hmap::Cmap::TERRAIN,
                  /* hillshading */ true);
  z_mise.to_png("docs/images/hydraulic_comparison_mise.png",
                hmap::Cmap::TERRAIN,
                /* hillshading */ true);

  hmap::export_banner_png("docs/images/hydraulic_comparison_terrains.png",
                          {z0, z_particle, z_stream, z_mise},
                          hmap::Cmap::TERRAIN,
                          /* hillshading */ true);

  hmap::export_banner_png("docs/images/hydraulic_comparison_erosion_maps.png",
                          {diff_particle, diff_stream, diff_mise},
                          hmap::Cmap::INFERNO,
                          /* hillshading */ false);

  std::cout << "Hydraulic comparison images successfully generated.\n";
  return 0;
}
