#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);

  // random seed points, elevation sampled from the heightmap
  const int   npoints = 5;
  const float radius = 0.2f;

  hmap::Cloud cloud(npoints, seed);
  cloud.set_values_from_array(z0, hmap::ArraySampling::AS_MEAN, 0.02f);

  const std::vector<float> xp = cloud.get_x();
  const std::vector<float> yp = cloud.get_y();
  const std::vector<float> zp = cloud.get_values();

  // terrain-aware parameters
  const float               slope_weight = 1.f;
  const float               uphill_weight = 2.f;
  const float               max_height_above = -1.f;
  const float               terrain_follow_plateau = 0.f;
  const float               terrain_follow_relief = 0.5f;
  const float               strength = 1.f;
  const hmap::RadialProfile radial_profile = hmap::RadialProfile::RP_SMOOTHSTEP;
  const float               profile_param = 1.f;

  // isotropic reference
  hmap::Array z_iso = z0;
  hmap::flatten_around_points(z_iso, xp, yp, zp, radius);

  // terrain-aware, exact plateau
  hmap::Array z_plateau = z0;
  hmap::Array mask_plateau;
  hmap::flatten_around_points_terrain_aware(z_plateau,
                                            xp,
                                            yp,
                                            zp,
                                            2.f * radius, // force a bit...
                                            slope_weight,
                                            uphill_weight,
                                            max_height_above,
                                            terrain_follow_plateau,
                                            strength,
                                            radial_profile,
                                            profile_param,
                                            &mask_plateau);

  // terrain-aware, keeping part of the relief
  hmap::Array z_follow = z0;
  hmap::Array mask_follow;
  hmap::flatten_around_points_terrain_aware(z_follow,
                                            xp,
                                            yp,
                                            zp,
                                            2.f * radius,
                                            slope_weight,
                                            uphill_weight,
                                            max_height_above,
                                            terrain_follow_relief,
                                            strength,
                                            radial_profile,
                                            profile_param,
                                            &mask_follow);

  hmap::export_banner_png("ex_flatten_around_points_terrain_aware.png",
                          {z0, z_iso, z_plateau, z_follow, mask_follow},
                          hmap::Cmap::TERRAIN,
                          true);
}
