#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};

  // surface defined by 4 corner elevations: c00, c10, c01, c11
  hmap::Array z1 = hmap::quad_surface(shape, 0.f, 1.f, 0.5f, 0.f);

  // another surface configuration
  hmap::Array z2 = hmap::quad_surface(shape, 1.f, 0.2f, 0.3f, 0.8f);

  hmap::export_banner_png("ex_quad_surface.png",
                          {z1, z2},
                          hmap::Cmap::INFERNO,
                          false);
}
