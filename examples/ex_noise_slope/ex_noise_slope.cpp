#include "highmap.hpp"

int main(void)
{
  hmap::NoiseType noise_type = hmap::NoiseType::SIMPLEX2;
  glm::ivec2      shape = {512, 512};
  glm::vec2       kw = {4.f, 4.f};
  std::uint32_t   seed = 1;
  int             octaves = 8;
  float           weight = 0.7f;
  float           persistence = 0.5f;
  float           lacunarity = 2.f;
  float           gradient_scale = 0.2f;
  float           slope_power = 2.f;
  float           octave_blend = 1.f;
  glm::vec4       bbox = {0.f, 1.f, 0.f, 1.f};
  float           angle = 30.f;
  float           stretch = 1.5f;
  float           stretch_decay = 0.5f;

  hmap::Array z = hmap::noise_slope(noise_type,
                                    shape,
                                    kw,
                                    seed,
                                    octaves,
                                    weight,
                                    persistence,
                                    lacunarity,
                                    gradient_scale,
                                    slope_power,
                                    octave_blend,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    bbox,
                                    angle,
                                    stretch,
                                    stretch_decay);

  z.to_png("ex_noise_slope.png", hmap::Cmap::TERRAIN, true);

  return 0;
}
