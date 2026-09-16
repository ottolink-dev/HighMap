#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 512};
  int        seed = 0;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {4.f, 4.f},
                                  seed);
  z.to_png("ex_generate_network_alpha_model0.png", hmap::Cmap::TERRAIN, true);

  // "cities"
  int         npoints = 4;
  int         ndummy = 50 * 50;
  float       alpha = 0.8f;
  glm::vec4   bbox = {1.f, 2.f, -0.5f, 0.5f};
  hmap::Cloud cloud = hmap::Cloud(npoints, seed, {1.2f, 1.8f, -0.4f, 0.4f});

  hmap::Graph network = hmap::generate_network_alpha_model(cloud.get_x(),
                                                           cloud.get_y(),
                                                           cloud.get_values(),
                                                           bbox,
                                                           z,
                                                           seed,
                                                           alpha,
                                                           ndummy);

  network.to_csv("nodes.csv", "adj.csv");

  // network.set_values(1.f);

  // graph view, no bounding box, color by edge weight
  network.to_png("ex_generate_network_alpha_model1.png");

  // graph project onto an array, with bounding box, color by point
  // values
  hmap::Array g = hmap::Array(shape);
  network.to_array(g, bbox);
  g.to_png("ex_generate_network_alpha_model2.png", hmap::Cmap::JET);
}
