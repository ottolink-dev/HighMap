#include "highmap.hpp"

int main(void)
{
  const glm::ivec2 shape = {256, 256};

  // 1. Define mountain ridge paths (A)

  hmap::Path ridge1;
  ridge1.add_point({0.150f, 0.200f, 1.f});
  ridge1.add_point({0.256f, 0.280f, 1.f});
  ridge1.add_point({0.360f, 0.220f, 1.f});
  ridge1 = hmap::catmullrom(ridge1);

  hmap::Path ridge2;
  ridge2.add_point({0.256f, 0.280f, 1.f});
  ridge2.add_point({0.280f, 0.380f, 1.f});
  ridge2 = hmap::catmullrom(ridge2);

  hmap::Array mountains(shape, 0.f);
  ridge1.to_array(mountains);
  ridge2.to_array(mountains);

  // 2. Define an optional coastline path (C)

  hmap::Path coast;
  int        n_coast_pts = 16;

  for (int i = 0; i <= n_coast_pts; ++i)
  {
    float angle = 2.f * M_PI * static_cast<float>(i) /
                  static_cast<float>(n_coast_pts);
    float r = 0.20f + 0.02f * std::sin(3.f * angle);

    coast.add_point(
        {0.256f + r * std::cos(angle), 0.256f + r * std::sin(angle), 1.f});
  }

  coast = hmap::catmullrom(coast);

  hmap::Array coastline(shape, 0.f);
  coast.to_array(coastline);

  // 3. Compute elevation using Rvachev R-function blending with different
  // exponents

  hmap::Array z_p05 = hmap::elevation_from_distance_fields(mountains,
                                                           nullptr,
                                                           &coastline,
                                                           0.5f);

  hmap::Array z_p10 = hmap::elevation_from_distance_fields(mountains,
                                                           nullptr,
                                                           &coastline,
                                                           1.f);

  hmap::Array z_p20 = hmap::elevation_from_distance_fields(mountains,
                                                           nullptr,
                                                           &coastline,
                                                           2.0f);

  // 4. Add procedural noise to distance fields

  hmap::Array noise = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {6.f, 6.f},
                                      42);

  hmap::Array z_noise = hmap::elevation_from_distance_fields(mountains,
                                                             nullptr,
                                                             &coastline,
                                                             1.5f,
                                                             1.f,
                                                             -1.f,
                                                             0.f,
                                                             &noise,
                                                             0.3f);

  hmap::export_banner_png("ex_elevation_from_distance_fields.png",
                          {z_p05, z_p10, z_p20, z_noise},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
