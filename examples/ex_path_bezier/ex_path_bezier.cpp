#include "highmap.hpp"

int main(void)
{
  glm::ivec2    shape = {256, 256};
  std::uint32_t seed = 1;

  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path path = hmap::Path(5, seed, hmap::adjust(bbox, -0.1f));
  path.reorder_nns();

  // --- Open path

  path.set_closed(false);
  auto z_sdf_open = hmap::path_sdf_to_array(path, shape, bbox);
  auto z_bezier_open = hmap::sdf_2d_polyline_bezier(path, shape, bbox);

  // --- Closed path

  path.set_closed(true);
  auto z_sdf_closed = hmap::path_sdf_to_array(path, shape, bbox);
  auto z_bezier_closed = hmap::sdf_2d_polyline_bezier(path, shape, bbox);

  hmap::export_banner_png(
      "ex_path_bezier.png",
      {z_sdf_open, z_bezier_open, z_sdf_closed, z_bezier_closed},
      hmap::Cmap::INFERNO);
}
