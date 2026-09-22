#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  shape = {512, 512};
  glm::vec2 res = {2.f, 2.f};
  int       seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);

  auto z1 = z;
  {
    int   ir = 16;
    float gamma = 2.f;
    hmap::recast_peak(z1, ir, gamma);
  }

  auto z2 = z;
  {
    hmap::Array noise = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                        shape,
                                        res,
                                        seed++);
    hmap::remap(noise, -0.1f, 0.1f);

    float zcut = 0.5f;
    float gamma = 4.f;
    hmap::recast_canyon(z2, zcut, gamma, &noise);
  }

  auto z3 = z;
  {
    float talus = 2.f / shape.x;
    int   ir = 4;
    float amplitude = 0.2f;
    float kw = 16.f;
    float gamma = 0.5f;
    hmap::recast_rocky_slopes(z3, talus, ir, amplitude, seed, kw, gamma);
  }

  auto z4 = z;
  auto z5 = z;
  {
    float talus = 1.f / shape.x;
    int   ir = 32;
    float amplitude = 0.05f;
    float gain = 1.5f;
    hmap::recast_cliff(z4, talus, ir, amplitude, gain);
    hmap::remap(z4);

    float angle = 0.f;
    hmap::recast_cliff_directional(z5, talus, ir, amplitude, angle, gain);
    hmap::remap(z5);
  }

  auto z6 = z;
  auto z7 = z;
  {
    float k = 0.01f;

    hmap::recast_billow(z6, 0.5f, k);
    hmap::remap(z6);

    hmap::recast_sag(z7, 0.5f, k);
    hmap::remap(z7);
  }

  auto z8 = z;
  {
    z8.to_png("out0.png", hmap::Cmap::GRAY, false, CV_16U);
    hmap::recast_escarpment(z8, 16, 0.1f, 1.f, false, true);

    z8.to_png("out1.png", hmap::Cmap::GRAY, false, CV_16U);
  }

  hmap::export_banner_png("ex_recast.png",
                          {z, z1, z2, z3, z4, z5, z6, z7, z8},
                          hmap::Cmap::TERRAIN,
                          true);
}
