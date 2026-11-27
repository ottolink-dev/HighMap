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

  hmap::Array upper_aligned = upper;
  hmap::align_terrains(lower,
                       frame_lower,
                       upper_aligned,
                       frame_upper,
                       hmap::TerrainAlignMode::TAM_PLANE);

  hmap::Array z_aligned_dist = lower;
  hmap::blend_terrains(z_aligned_dist,
                       frame_lower,
                       upper_aligned,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_DISTANCE,
                       0.2f);

  hmap::Array z_height = lower;
  hmap::blend_terrains(z_height,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_HEIGHT,
                       0.2f);

  hmap::Array z_slope = lower;
  hmap::blend_terrains(z_slope,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_SLOPE,
                       0.2f);

  hmap::Array z_seamless = lower;
  hmap::blend_terrains(z_seamless,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_SEAMLESS,
                       0.2f);

  hmap::Array z_multiscale = lower;
  hmap::blend_terrains(z_multiscale,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_MULTISCALE,
                       0.25f);

  hmap::Array z_poisson = lower;
  hmap::blend_terrains(z_poisson,
                       frame_lower,
                       upper,
                       frame_upper,
                       hmap::TerrainBlendMode::TBM_POISSON,
                       0.15f);

  hmap::export_banner_png("ex_blend_terrains0.png",
                          {lower, z_max, z_max_smooth, z_dist, z_aligned_dist},
                          hmap::Cmap::TERRAIN);

  hmap::export_banner_png(
      "ex_blend_terrains1.png",
      {z_height, z_slope, z_seamless, z_multiscale, z_poisson},
      hmap::Cmap::TERRAIN);
}
