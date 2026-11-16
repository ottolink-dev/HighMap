#include "highmap.hpp"
#include "highmap/internal/timer.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {6.f, 6.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);

  float additional_depth = 0.05f;

  hmap::Array water_depth0 = hmap::flooding_uniform_level(z, 0.3f);
  // hmap::Array water_depth0 = hmap::flooding_lake_system(z);

  hmap::Array water_depth1 = hmap::water_depth_increase(water_depth0,
                                                        z,
                                                        additional_depth);

  hmap::export_banner_png("ex_water_depth_increase.png",
                          {z, z + water_depth0, z + water_depth1},
                          hmap::Cmap::TERRAIN,
                          true);
}
