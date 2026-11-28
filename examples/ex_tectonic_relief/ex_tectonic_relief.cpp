#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};

  int   plate_count = 32;
  uint  seed = 1;
  float continental_ratio = 0.2f;
  float smoothing_coherence = 0.05f;
  float roughness = 1.f;
  float boundary_width = 0.08f;
  float relief_amp = 1.0f;
  float subduction_offset = 0.8f;
  float fault_amp = 0.4f;
  float ripple_amp = 0.1f;
  float ripple_kw = 96.0f;
  bool  underwater_interactions = false;

  hmap::Array ranges;
  hmap::Array faults;
  hmap::Array land;

  glm::vec2   kw = {8.f, 8.f};
  hmap::Array noise = 1.5f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                             shape,
                                             kw,
                                             seed + 1);

  hmap::Array z = hmap::tectonic_relief(shape,
                                        plate_count,
                                        seed,
                                        continental_ratio,
                                        smoothing_coherence,
                                        roughness,
                                        boundary_width,
                                        relief_amp,
                                        subduction_offset,
                                        fault_amp,
                                        ripple_amp,
                                        ripple_kw,
                                        underwater_interactions,
                                        nullptr,
                                        nullptr,
                                        &noise,
                                        &noise,
                                        &ranges,
                                        &faults,
                                        &land);

  z.dump();

  hmap::export_banner_png("ex_tectonic_relief.png",
                          {z, ranges, faults, land},
                          hmap::Cmap::TERRAIN);
}
