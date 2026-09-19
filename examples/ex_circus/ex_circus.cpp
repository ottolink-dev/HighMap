#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 512};

  // --- 1. Default / baseline circus envelope (exit towards east)

  float radius = 0.35f;
  float angle = 0.f;
  float exit_width = 0.35f;
  float exit_depth = 0.0f;
  float center_depth = 0.1f;
  float ridge_height = 1.0f;
  float ridge_width = 0.15f;
  float outer_falloff = 2.0f;

  hmap::Array z1 = hmap::circus(shape,
                                radius,
                                angle,
                                exit_width,
                                exit_depth,
                                center_depth,
                                ridge_height,
                                ridge_width,
                                outer_falloff);

  // --- 2. Wide exit angle, rotated orientation, deeper depression

  float angle_var = 45.f;
  float exit_width_var = 0.55f;
  float center_depth_var = 0.02f;
  float ridge_height_var = 0.85f;
  float ridge_width_var = 0.25f;

  hmap::Array z2 = hmap::circus(shape,
                                radius,
                                angle_var,
                                exit_width_var,
                                exit_depth,
                                center_depth_var,
                                ridge_height_var,
                                ridge_width_var,
                                outer_falloff);

  // --- 3. Circus modulatedwith procedural noise

  glm::vec2     kw = {4.f, 4.f};
  std::uint32_t seed = 42;

  auto noise_r = 0.15f *
                 hmap::noise_fbm(hmap::NoiseType::SIMPLEX2S, shape, kw, seed);

  auto terrain_noise = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2S,
                                       shape,
                                       kw * 2.f,
                                       seed + 1);
  hmap::remap(terrain_noise, 0.f, 0.05f);

  hmap::Array z3 = hmap::circus(shape,
                                radius,
                                -30.f,
                                exit_width,
                                exit_depth,
                                center_depth,
                                ridge_height,
                                ridge_width,
                                outer_falloff,
                                &noise_r);

  // --- Export comparison banner

  hmap::export_banner_png("ex_circus.png",
                          {z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
