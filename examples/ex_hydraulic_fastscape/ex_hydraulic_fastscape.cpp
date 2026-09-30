#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  // shape = {1024, 1024};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 0;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z0);
  hmap::Array z = z0;
  hmap::Array erosion_map(shape, 0.f);
  hmap::Array flow_map(shape, 0.f);

  // Apply FastScape landscape evolution
  hmap::hydraulic_fastscape(z,
                            50,      // iterations
                            1e-2f,   // dt
                            1.f,     // k_erosion
                            0.5f,    // m_exp
                            1.f,     // n_exp
                            1e-3f,   // k_diff
                            0.1f,   // uplift_rate
                            true,    // multiple_flow
                            1.f,     // flow_partition_exp
                            1e-3f,   // tolerance
                            nullptr, // bedrock
                            nullptr, // moisture
                            &erosion_map,
                            &flow_map);

  z.dump();
  
  hmap::remap(erosion_map);

  hmap::export_banner_png("ex_hydraulic_fastscape.png",
                          {z0, z, erosion_map},
                          hmap::Cmap::TERRAIN,
                          true,
                          true);
  return 0;
}
