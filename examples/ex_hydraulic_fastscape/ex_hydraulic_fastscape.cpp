#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  // shape = {1024, 1024};
  // shape = {512, 512};
  glm::vec2 res = {2.f, 2.f};
  int       seed = 0;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z0);
  hmap::Array z = z0;
  hmap::Array erosion_map(shape, 0.f);
  hmap::Array flow_map(shape, 0.f);

  // Shared erosion parameters
  float dt = 1e-2f;
  float k_erosion = 1.f;
  float m_exp = 0.5f;
  float n_exp = 0.8f;
  float k_diff = 1e-3f;
  float uplift_rate = 0.1f;
  bool  multiple_flow = true;
  float flow_partition_exp = 1.f;
  float tolerance = 1e-3f;

  // Apply FastScape landscape evolution (single-scale)
  int iterations = 16;
  hmap::hydraulic_fastscape(z,
                            iterations,
                            dt,
                            k_erosion,
                            m_exp,
                            n_exp,
                            k_diff,
                            uplift_rate,
                            multiple_flow,
                            flow_partition_exp,
                            tolerance,
                            nullptr, // bedrock
                            nullptr, // moisture
                            &erosion_map,
                            &flow_map);

  // Apply Multiscale FastScape landscape evolution
  hmap::Array z_multi = z0;
  hmap::Array erosion_map_multi(shape, 0.f);
  hmap::Array flow_map_multi(shape, 0.f);

  std::vector<int> steps_per_level = {16, 8, 4};
  float            mix = 1.f;

  hmap::hydraulic_fastscape_multiscale(z_multi,
                                       steps_per_level,
                                       dt,
                                       k_erosion,
                                       m_exp,
                                       n_exp,
                                       k_diff,
                                       uplift_rate,
                                       multiple_flow,
                                       flow_partition_exp,
                                       tolerance,
                                       nullptr, // bedrock
                                       nullptr, // moisture
                                       &erosion_map_multi,
                                       &flow_map_multi,
                                       mix);

  hmap::remap(erosion_map);
  hmap::remap(erosion_map_multi);

  z_multi.dump();
  lerp(z0, z_multi, 0.5f).dump("out1.png");
  
  hmap::export_banner_png("ex_hydraulic_fastscape.png",
                          {z0, z, erosion_map, z_multi, erosion_map_multi},
                          hmap::Cmap::TERRAIN,
                          true,
                          true);
  return 0;
}
