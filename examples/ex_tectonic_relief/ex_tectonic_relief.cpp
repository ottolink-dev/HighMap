#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        plate_count = 32;
  uint       seed = 0;
  float      continental_ratio = 0.4f;
  float      roughness = 0.1f;
  float      boundary_width = 0.08f;
  float      relief_amp = 1.0f;
  float      subduction_offset = 0.8f;
  float      fault_amp = 0.15f;

  hmap::Array ranges;
  hmap::Array faults;
  hmap::Array land;

  hmap::Array z = hmap::tectonic_relief(shape,
                                        plate_count,
                                        seed,
                                        continental_ratio,
                                        roughness,
                                        boundary_width,
                                        relief_amp,
                                        subduction_offset,
                                        fault_amp,
                                        nullptr,
                                        nullptr,
                                        nullptr,
                                        nullptr,
                                        &ranges,
                                        &faults,
                                        &land);

  hmap::export_banner_png("ex_tectonic_relief.png",
                          {z, ranges, faults, land},
                          hmap::Cmap::TERRAIN);
}
