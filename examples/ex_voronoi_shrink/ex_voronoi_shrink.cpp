#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  const glm::ivec2 shape = {512, 512};
  const glm::vec2  kw_base = {2.f, 2.f};
  const int        seed = 42;

  // Base terrain heightmap (continent-scale terrain)
  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw_base,
                                   seed);
  hmap::remap(z0, 0.f, 1.f);

  const float     kw = 4.f;
  const glm::vec2 jitter = {0.6f, 0.6f};
  const float     fill_ocean = 0.f;

  // Mild shrinking (narrow rifts)
  auto z1 = hmap::gpu::voronoi_shrink(z0, kw, 0.85f, fill_ocean, seed, jitter);

  // Moderate shrinking (clearly separated continents with matching coastlines)
  auto z2 = hmap::gpu::voronoi_shrink(z0, kw, 0.70f, fill_ocean, seed, jitter);

  // Strong shrinking with domain noise displacement on rifts
  hmap::Array noise_x = 0.2f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                               shape,
                                               {4.f, 4.f},
                                               seed + 1);
  hmap::Array noise_y = 0.2f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                               shape,
                                               {4.f, 4.f},
                                               seed + 2);

  // Point-based Voronoi shrink with custom seed centers
  std::vector<glm::vec2> custom_centers = {{0.25f, 0.3f},
                                           {0.75f, 0.25f},
                                           {0.3f, 0.75f},
                                           {0.7f, 0.7f},
                                           {0.5f, 0.45f}};

  auto z3 = hmap::gpu::voronoi_shrink(z0,
                                      custom_centers,
                                      0.65f,
                                      fill_ocean,
                                      nullptr,
                                      &noise_x,
                                      &noise_y);

  hmap::export_banner_png("ex_voronoi_shrink.png",
                          {z0, z1, z2, z3},
                          hmap::Cmap::TERRAIN,
                          true);
}
