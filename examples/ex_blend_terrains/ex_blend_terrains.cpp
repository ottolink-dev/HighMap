#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 1;

  hmap::Array lower = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {2.f, 2.f},
                                      seed);
  hmap::Array upper = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      {128, 128},
                                      {4.f, 4.f},
                                      seed + 1);

  hmap::remap(lower, 0.f, 0.5f);
  hmap::remap(upper, 0.2f, 1.0f);

  hmap::CoordFrame frame_lower({0.f, 0.f}, {1.f, 1.f}, 0.f);
  hmap::CoordFrame frame_upper({0.25f, 0.2f}, {0.55f, 0.55f}, 15.f);

  hmap::Array z_add = lower;
  hmap::blend_terrains(z_add,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_ADD);

  hmap::Array z_max = lower;
  hmap::blend_terrains(z_max,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_MAX);

  hmap::Array z_max_smooth = lower;
  hmap::blend_terrains(z_max_smooth,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_MIN,
                       0.25f);

  hmap::Array z_dist = lower;
  hmap::blend_terrains(z_dist,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_DISTANCE,
                       0.2f);

  hmap::export_banner_png("ex_blend_terrains0.png",
                          {lower, z_add, z_max, z_max_smooth, z_dist},
                          hmap::Cmap::TERRAIN);
}
