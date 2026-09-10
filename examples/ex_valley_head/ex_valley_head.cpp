#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 512};

  // --- Noise for cross-section and width modulation

  glm::vec2     kw = {6.f, 6.f};
  std::uint32_t seed = 42;

  // displacement noise on local cross-section (meander)
  auto noise_offset = 0.3f * hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                             shape,
                                             kw,
                                             ++seed);

  // radial / width noise
  auto noise_r = 0.1f * hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                        shape,
                                        kw * 2.f,
                                        ++seed);

  // --- Clean valley head incision

  float     angle = 180.f; // along -y
  float     max_depth = 0.25f;
  float     min_width = 0.02f;
  float     max_width = 0.5f;
  float     depth_length = 0.45f;
  float     width_length = 0.45f;
  float     dev_power = 2.f;
  float     profile_power = 1.3f;
  glm::vec2 center = {0.5f, 1.f};

  auto delta_clean = hmap::valley_head(shape,
                                       angle,
                                       max_depth,
                                       min_width,
                                       max_width,
                                       depth_length,
                                       width_length,
                                       dev_power,
                                       profile_power,
                                       nullptr,
                                       nullptr,
                                       center);

  // --- Noisy valley head incision

  auto delta_noisy = hmap::valley_head(shape,
                                       angle,
                                       max_depth,
                                       min_width,
                                       max_width,
                                       depth_length,
                                       width_length,
                                       dev_power,
                                       profile_power,
                                       &noise_offset,
                                       &noise_r,
                                       center);

  // --- Export

  hmap::export_banner_png("ex_valley_head.png",
                          {delta_clean, delta_noisy},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
