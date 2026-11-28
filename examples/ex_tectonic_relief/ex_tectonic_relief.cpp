#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};

  hmap::TectonicReliefParams params;
  params.plate_count = 32;
  params.seed = 0;
  params.boundary_width = 0.08f;
  params.subduction_offset = 0.1f;

  hmap::Array ranges;
  hmap::Array faults;
  hmap::Array z = hmap::tectonic_relief(shape,
                                        params,
                                        nullptr,
                                        nullptr,
                                        &ranges,
                                        &faults);

  auto zc = z;
  hmap::clamp_min(zc, 0.f);

  hmap::export_banner_png("ex_tectonic_relief.png",
                          {z, ranges, faults, zc},
                          hmap::Cmap::TERRAIN);
}
