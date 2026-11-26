#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 3;

  hmap::BBox2D bbox = hmap::unit_square_bbox();
  hmap::Path   path = hmap::Path(4, seed, bbox);
  path.reorder_nns();
  path = hmap::bezier(path);

  // input array, scattered non-zero values
  hmap::Array z0 = hmap::Array(shape);
  path.to_array(z0, bbox);

  hmap::Array z1 = reverse_midpoint(z0, seed++);
  hmap::median_3x3(z1);

  hmap::export_banner_png("ex_reverse_midpoint.png", {z0, z1}, hmap::Cmap::JET);
}
