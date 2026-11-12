#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {1024, 1024};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);

  // random cloud, with values sampled from the heightmap (min in a radius)
  const int   npoints = 5;
  const float sampling_radius = 0.1f;
  const float flatten_radius = 0.2f;
  hmap::Cloud cloud(npoints, seed);

  cloud.set_values_from_array(z0,
                              hmap::ArraySampling::AS_MEAN,
                              sampling_radius);

  const std::vector<float> xp = cloud.get_x();
  const std::vector<float> yp = cloud.get_y();
  const std::vector<float> zp = cloud.get_values();

  hmap::Array z1 = z0;
  hmap::flatten_around_points(z1, xp, yp, zp, flatten_radius);

  hmap::export_banner_png("ex_flatten_around_points.png",
                          {z0, z1},
                          hmap::Cmap::TERRAIN,
                          true);
}
