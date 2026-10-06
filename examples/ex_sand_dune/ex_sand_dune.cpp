#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {3.f, 3.f};
  int        seed = 42;

  // base terrain with rolling bedrock hills and initial sand layer
  hmap::Array bedrock = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(bedrock, 0.1f, 0.6f);

  hmap::Array z0 = bedrock + 0.01f; // initial sand mantle
  hmap::Array z_dunes = z0;

  hmap::SandDuneParams params;
  params.wind_angle = 0.f;            // wind blowing predominantly eastward
  params.wind_speed = 1.2f;
  params.crest_speedup = 0.6f;         // crest acceleration
  params.hop_length = 6.0f;
  params.elevation_hop_factor = 3.0f;  // elevation-aware ballistic hops
  params.jitter_angle = 0.15f;
  params.p_sand = 0.65f;
  params.p_bare = 0.35f;
  params.shadow_talus = 1.5f/ shape.x;        // shadow height difference threshold between 2 cells
  params.talus = 2.f / shape.x;               // angle of repose talus (height difference between 2 cells)
  params.slab_height = 0.001f;
  params.collapse_rate = 0.5f;
  params.periodic_boundary = true;     // toroidal domain

  hmap::Array dep_map(shape, 0.f);
  hmap::Array ero_map(shape, 0.f);

  int iterations = 1024;
  
  hmap::log::info("Simulating sand dunes and aeolian deposition...");
  hmap::gpu::sand_dune(z_dunes,
                       500000,
                       params,
                       &bedrock,
                       nullptr,
                       &dep_map,
                       &ero_map,
                       iterations);

  hmap::log::info("z0 elevation range: [{:.3f}, {:.3f}]", z0.min(), z0.max());
  hmap::log::info("z_dunes elevation range: [{:.3f}, {:.3f}]", z_dunes.min(), z_dunes.max());

  hmap::export_banner_png("ex_sand_dune.png",
                          {z0, z_dunes, dep_map},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
