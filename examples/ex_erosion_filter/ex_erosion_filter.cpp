#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {2.f, 2.f};
  int        seed = 1;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                   shape,
                                   kw,
                                   seed,
                                   /* octaves */ 2,
                                   /* weight */ 1.f,
                                   /* persistence */ 0.5f,
                                   /* lacunarity */ 2.f);
  hmap::remap(z0, 0.f, 1.f);

  // Apply Runevision Advanced Terrain Erosion Filter
  hmap::Array z_eroded = z0;
  hmap::Array ridge_map(shape);

  hmap::gpu::erosion_filter(z_eroded,
                            /* scale */ 0.15f,
                            /* strength */ 0.22f,
                            /* gully_weight */ 2.f,
                            /* detail */ 1.5f,
                            /* rounding */ {0.5f, 0.5f, 1.f, 1.f},
                            /* onset */ {1.25f, 0.5f, 1.f, 0.5f},
                            /* assumed_slope */ {4.f, 1.f},
                            /* cell_scale */ 0.8f,
                            /* octaves */ 8,
                            /* gain */ 0.7f,
                            /* lacunarity */ 2.0f,
                            /* normalization */ 0.4f,
                            /* curvature_scaling */ 10.f,
                            /* seed */ 0,
                            /* p_fade_target */ nullptr,
                            /* p_ridge_map */ &ridge_map);

  hmap::remap(z_eroded, 0.f, 1.f);

  hmap::Array ridge_vis = ridge_map;
  hmap::remap(ridge_vis, 0.f, 1.f);

  hmap::export_banner_png("ex_erosion_filter.png",
                          {z0, z_eroded, ridge_vis},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
