#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {4.f, 4.f};
  int        seed = 0;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z0);

  float riverbed_talus = 0.01f / shape.x;

  auto z1 = hmap::flow_fixing(z0, riverbed_talus);
  auto z2 = hmap::flow_fixing_drainage_basin(z0,
                                             hmap::FlowDirectionMethod::FDM_D8,
                                             riverbed_talus,
                                             50,
                                             true);
  auto z3 = hmap::flow_fixing_mst(z0,
                                  riverbed_talus,
                                  0.99f, // elevation_ratio
                                  2.f,   // distance_exponent
                                  100.f, // upward_penalization
                                  0.5f,  // valley_affinity
                                  8,     // prefilter_ir
                                  1e-4f, // minimum_depth
                                  true,  // carve_riverbed
                                  8.f,   // merging_distance
                                  hmap::RadialProfile::RP_SMOOTHSTEP_UPPER,
                                  2.f,     // radial_profile_parameter
                                  nullptr, // p_noise_r
                                  0,       // fractalize_iterations
                                  0.1f,    // fractalize_sigma
                                  4,       // decimate_target_points
                                  seed);   // fractalize_seed

  z3.dump();

  hmap::export_banner_png("ex_flow_fixing.png",
                          {z0, z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
