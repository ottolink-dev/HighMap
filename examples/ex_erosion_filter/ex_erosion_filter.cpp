#include <vector>

#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {2.f, 2.f};
  int        seed = 1;

  // Generate input base terrain
  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                   shape,
                                   kw,
                                   seed,
                                   /* octaves */ 1,
                                   /* weight */ 0.7f,
                                   /* persistence */ 0.5f,
                                   /* lacunarity */ 2.f);
  hmap::remap(z0, 0.f, 1.f);

  std::vector<hmap::Array> banner = {z0};

  // --- Varying strength (erosion depth)

  for (float strength : {0.05f, 0.25f, 0.45f})
  {
    hmap::Array z = z0;
    hmap::gpu::erosion_filter(z,
                              /* scale */ 0.15f,
                              strength,
                              /* gully_weight */ 0.5f,
                              /* detail */ 1.f,
                              /* rounding */ {0.5f, 0.5f, 1.f, 1.f},
                              /* onset */ {1.25f, 0.5f, 1.f, 0.5f},
                              /* assumed_slope */ {4.f, 1.f},
                              /* cell_scale */ 0.8f,
                              /* octaves */ 8,
                              /* gain */ 0.7f,
                              /* lacunarity */ 2.0f,
                              /* normalization */ 0.4f,
                              /* curvature_scaling */ 10.f,
                              /* seed */ 0);
    hmap::remap(z, 0.f, 1.f);
    banner.push_back(z);

    if (strength == 0.25f) z.dump();
  }

  // --- Varying scale (feature size)

  for (float scale : {0.08f, 0.30f})
  {
    hmap::Array z = z0;
    hmap::gpu::erosion_filter(z,
                              scale,
                              /* strength */ 0.25f,
                              /* gully_weight */ 0.5f,
                              /* detail */ 1.f,
                              /* rounding */ {0.5f, 0.5f, 1.f, 1.f},
                              /* onset */ {1.25f, 0.5f, 1.f, 0.5f},
                              /* assumed_slope */ {4.f, 1.f},
                              /* cell_scale */ 0.8f,
                              /* octaves */ 8,
                              /* gain */ 0.7f,
                              /* lacunarity */ 2.0f,
                              /* normalization */ 0.4f,
                              /* curvature_scaling */ 10.f,
                              /* seed */ 0);
    hmap::remap(z, 0.f, 1.f);
    banner.push_back(z);
  }

  // --- Varying gully weight (gully carving vs diffuse erosion)

  for (float gully_weight : {0.5f, 1.f, 3.f})
  {
    hmap::Array z = z0;
    hmap::gpu::erosion_filter(z,
                              /* scale */ 0.15f,
                              /* strength */ 0.25f,
                              gully_weight,
                              /* detail */ 1.f,
                              /* rounding */ {0.5f, 0.5f, 1.f, 1.f},
                              /* onset */ {1.25f, 0.5f, 1.f, 0.5f},
                              /* assumed_slope */ {4.f, 1.f},
                              /* cell_scale */ 0.8f,
                              /* octaves */ 8,
                              /* gain */ 0.7f,
                              /* lacunarity */ 2.0f,
                              /* normalization */ 0.4f,
                              /* curvature_scaling */ 10.f,
                              /* seed */ 0);
    hmap::remap(z, 0.f, 1.f);
    banner.push_back(z);
  }

  // --- Ridge map output

  hmap::Array z_eroded = z0;
  hmap::Array ridge_map(shape);
  hmap::gpu::erosion_filter(z_eroded,
                            /* scale */ 0.15f,
                            /* strength */ 0.25f,
                            /* gully_weight */ 0.5f,
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
  hmap::remap(ridge_map, 0.f, 1.f);
  banner.push_back(ridge_map);

  // Export comparison banner
  hmap::export_banner_png("ex_erosion_filter.png",
                          banner,
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
