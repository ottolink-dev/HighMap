#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};

  // 1. Define mountain ridge lines (values > 0)
  hmap::Path ridge1;
  ridge1.push_back({0.150f, 0.250f, 1.0f});
  ridge1.push_back({0.260f, 0.230f, 1.0f});
  ridge1.push_back({0.360f, 0.220f, 1.0f});
  ridge1 = hmap::catmullrom(ridge1);

  hmap::Path ridge2;
  ridge2.push_back({0.256f, 0.280f, 1.0f});
  ridge2.push_back({0.280f, 0.380f, 1.0f});
  ridge2 = hmap::catmullrom(ridge2);

  hmap::Array mountains(shape, 0.0f);
  ridge1.to_array(mountains);
  ridge2.to_array(mountains);

  // 2. Define closed coastline path (z = 0)
  hmap::Path coast;
  int        n_coast_pts = 24;

  for (int i = 0; i <= n_coast_pts; ++i)
  {
    float angle = 2.0f * M_PI * static_cast<float>(i) /
                  static_cast<float>(n_coast_pts);
    float r = 0.22f + 0.03f * std::sin(3.0f * angle);

    coast.push_back(
        {0.256f + r * std::cos(angle), 0.256f + r * std::sin(angle), 1.0f});
  }
  coast = hmap::catmullrom(coast);

  hmap::Array coastline(shape, 0.0f);
  coast.to_array(coastline);

  // 3. Generate smooth continuous heightmap using GPU harmonic interpolation
  hmap::Array z_harmonic = hmap::elevation_from_sparse_constraints(mountains,
                                                                   &coastline,
                                                                   nullptr,
                                                                   1.0f,
                                                                   0.0f,
                                                                   -1.0f);

  // 4. Modulate with optional procedural noise
  hmap::Array noise = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {6.0f, 6.0f},
                                      42);

  hmap::Array z_detailed = hmap::elevation_from_sparse_constraints(mountains,
                                                                   &coastline,
                                                                   nullptr,
                                                                   1.0f,
                                                                   0.0f,
                                                                   -1.0f,
                                                                   500,
                                                                   1e-5f,
                                                                   &noise,
                                                                   0.25f);

  hmap::export_banner_png("ex_elevation_from_sparse_constraints.png",
                          {mountains, coastline, z_harmonic, z_detailed},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
