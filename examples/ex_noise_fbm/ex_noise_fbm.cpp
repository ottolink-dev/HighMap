#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  std::vector<hmap::NoiseType> noise_list = {hmap::NoiseType::PERLIN,
                                             hmap::NoiseType::PERLIN_BILLOW,
                                             hmap::NoiseType::PERLIN_HALF,
                                             hmap::NoiseType::SIMPLEX2,
                                             hmap::NoiseType::SIMPLEX2S,
                                             hmap::NoiseType::VALUE,
                                             hmap::NoiseType::VALUE_CUBIC,
                                             hmap::NoiseType::VALUE_DELAUNAY,
                                             hmap::NoiseType::VALUE_LINEAR,
                                             hmap::NoiseType::WORLEY,
                                             hmap::NoiseType::WORLEY_DOUBLE,
                                             hmap::NoiseType::WORLEY_VALUE};

  auto ctrl_array = hmap::noise(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(ctrl_array);

  // fbm
  {
    std::vector<hmap::Array> arrays;
    for (auto &noise_type : noise_list)
    {
      auto n = hmap::noise_fbm(noise_type,
                               shape,
                               kw,
                               seed,
                               8,
                               0.7f,
                               0.5f,
                               2.f,
                               &ctrl_array);
      hmap::remap(n);
      arrays.push_back(n);
    }
    hmap::export_banner_png("ex_noise_fbm0.png",
                            arrays,
                            hmap::Cmap::TERRAIN,
                            true);
  }

  // iq
  {
    std::vector<hmap::Array> arrays;
    for (auto &noise_type : noise_list)
    {
      auto n = hmap::noise_iq(noise_type, shape, kw, seed);
      hmap::remap(n);
      arrays.push_back(n);
    }
    hmap::export_banner_png("ex_noise_fbm1.png",
                            arrays,
                            hmap::Cmap::TERRAIN,
                            true);
  }

  // jordan
  {
    std::vector<hmap::Array> arrays;
    for (auto &noise_type : noise_list)
    {
      auto n = hmap::noise_jordan(noise_type, shape, kw, seed);
      hmap::remap(n);
      arrays.push_back(n);
    }
    hmap::export_banner_png("ex_noise_fbm2.png",
                            arrays,
                            hmap::Cmap::TERRAIN,
                            true);
  }

  // pingpong
  {
    std::vector<hmap::Array> arrays;
    for (auto &noise_type : noise_list)
    {
      auto n = hmap::noise_pingpong(noise_type, shape, kw, seed);
      hmap::remap(n);
      arrays.push_back(n);
    }
    hmap::export_banner_png("ex_noise_fbm3.png",
                            arrays,
                            hmap::Cmap::TERRAIN,
                            true);
  }

  // ridged
  {
    std::vector<hmap::Array> arrays;
    for (auto &noise_type : noise_list)
    {
      auto n = hmap::noise_ridged(noise_type, shape, kw, seed);
      hmap::remap(n);
      arrays.push_back(n);
    }
    hmap::export_banner_png("ex_noise_fbm4.png",
                            arrays,
                            hmap::Cmap::TERRAIN,
                            true);
  }

  // swiss
  {
    std::vector<hmap::Array> arrays;
    for (auto &noise_type : noise_list)
    {
      auto n = hmap::noise_swiss(noise_type, shape, kw, seed);
      hmap::remap(n);
      arrays.push_back(n);
    }
    hmap::export_banner_png("ex_noise_fbm5.png",
                            arrays,
                            hmap::Cmap::TERRAIN,
                            true);
  }

  // Parberry
  {
    auto n = hmap::noise_parberry(shape, kw, seed);
    hmap::remap(n);
    n.to_png("ex_noise_fbm6.png", hmap::Cmap::TERRAIN, true);
  }

  // stretched fbm (isotropic, stretched, stretched with octave decay)
  {
    auto n0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
    auto n1 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                              shape,
                              kw,
                              seed,
                              8,
                              0.7f,
                              0.5f,
                              2.f,
                              nullptr,
                              nullptr,
                              nullptr,
                              {0.f, 1.f, 0.f, 1.f},
                              30.f,
                              4.f,
                              0.f);
    auto n2 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                              shape,
                              kw,
                              seed,
                              8,
                              0.7f,
                              0.5f,
                              2.f,
                              nullptr,
                              nullptr,
                              nullptr,
                              {0.f, 1.f, 0.f, 1.f},
                              30.f,
                              4.f,
                              0.5f);

    hmap::remap(n0);
    hmap::remap(n1);
    hmap::remap(n2);

    hmap::export_banner_png("ex_noise_fbm7.png",
                            {n0, n1, n2},
                            hmap::Cmap::TERRAIN,
                            true);
  }
}
