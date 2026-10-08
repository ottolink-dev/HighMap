#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  glm::vec2  kw = {8.f, 8.f};
  int        seed = 1;
  glm::vec2  jitter = {1.f, 1.f};

  hmap::Array z1 = hmap::gpu::voronoi_edge_distance(shape, kw, seed, jitter);

  hmap::Array noise = 0.2f * hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                             shape,
                                             {4.f, 4.f},
                                             seed + 1);

  hmap::Array z2 = hmap::gpu::voronoi_edge_distance(shape,
                                                    kw,
                                                    seed,
                                                    jitter,
                                                    /* bias */ 0.f,
                                                    nullptr,
                                                    &noise);

  z2.infos();

  hmap::export_banner_png("ex_voronoi_edge_distance.png",
                          {z1, z2},
                          hmap::Cmap::INFERNO);
}
