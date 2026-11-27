#include "highmap.hpp"

int main(void)
{
  glm::ivec2    shape = {256, 256};
  std::uint32_t seed = 1;

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path path = hmap::Path(15, seed, hmap::BBox2D(0.1f, 0.9f, 0.1f, 0.9f));
  path.reorder_nns();

  auto z0 = path.to_array(shape, bbox);
  auto z1 = hmap::path_sdf_to_array(path, shape, bbox);

  path.set_closed(true);
  auto z2 = hmap::path_sdf_to_array(path, shape, bbox);

  hmap::export_banner_png("ex_path_sdf.png", {z0, z1, z2}, hmap::Cmap::INFERNO);
}
