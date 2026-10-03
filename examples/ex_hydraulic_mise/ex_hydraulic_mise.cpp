#include <iostream>

#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {3.f, 3.f};
  int        seed = 42;

  // --- Initial heightmap generation

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z);
  hmap::Array z0 = z;

  // --- Mise erosion solver configuration

  hmap::MiseParams params;
  params.strength = 0.6f;
  params.area_exp = 0.45f;
  params.talus = 3.2f;
  params.debris = 200.f;
  params.deposition_rate = 1.f;
  params.base_res = 128;
  params.seed = seed;

  hmap::Array sediment;
  hmap::Array flow;
  hmap::Array erosion_map;
  hmap::Array deposition_map;

  // run multiscale implicit stream-power erosion
  hmap::hydraulic_mise(z,
                       params,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       &sediment,
                       &flow,
                       &erosion_map,
                       &deposition_map);

  // --- Output and visualization

  std::cout << "original z range: [" << z0.min() << ", " << z0.max() << "]\n";
  std::cout << "eroded z range:   [" << z.min() << ", " << z.max() << "]\n";
  std::cout << "sediment max:     " << sediment.max() << "\n";
  std::cout << "erosion max:      " << erosion_map.max() << "\n";
  std::cout << "deposition max:   " << deposition_map.max() << "\n";

  // export before/after comparison banner
  hmap::export_banner_png("ex_hydraulic_mise.png",
                          {z0, z, sediment, flow},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
