#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  int        seed = 0;

  // --- generate mask

  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  hmap::Path path = hmap::Path(10, seed, hmap::BBox2D(0.1f, 0.9f, 0.1f, 0.9f));
  path.reorder_nns();
  auto mask = hmap::bezier(path).to_array(shape, bbox);
  hmap::clamp_min(mask, 0.f);

  // --- splat kernel

  glm::ivec2  kernel_shape = {64, 64};
  hmap::Array kernel = hmap::cubic_pulse(kernel_shape);

  auto z = hmap::gpu::sparse_max_convolution(mask, kernel);

  hmap::export_banner_png("ex_sparse_max_convolution.png",
                          {mask, z},
                          hmap::Cmap::JET);
}
