#include "highmap.hpp"

int main(void)
{
  hmap::log::warn(
      "hydraulic_stream_upscale_amplification is deprecated and will be "
      "removed at some point.");

  glm::ivec2 shape = {256, 256};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);
  auto z0 = z;

  float c_erosion = 0.05f;
  float talus_ref = 5.f / (float)shape.x;

  auto z1 = z;
  hmap::hydraulic_stream(z1, c_erosion, talus_ref);

  auto z2 = z;
  int  upscaling_levels = 2;
  hmap::hydraulic_stream_upscale_amplification(z2,
                                               c_erosion,
                                               talus_ref,
                                               upscaling_levels);

  auto z3 = z;
  upscaling_levels = 3;
  float persistence = 0.75f; // decrease ratio of erosion coeff at each level
  hmap::hydraulic_stream_upscale_amplification(z3,
                                               c_erosion,
                                               talus_ref,
                                               upscaling_levels,
                                               persistence);

  hmap::export_banner_png("ex_hydraulic_stream_upscale_amplification.png",
                          {z0, z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
