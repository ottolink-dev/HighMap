#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 1;

  hmap::BBox2D bbox = {-1.f, 2.f, 0.f, 5.f};

  // --- open path with (x, y) and values defined as vectors

  hmap::Path path = hmap::Path({-0.5f, 1.5f, 0.5f}, // x
                               {1.f, 2.f, 4.f},     // y
                               {0.f, 1.f, 2.f}      // values
  );

  hmap::Array z1 = hmap::Array(shape);
  path.to_array(z1, bbox);

  // --- generate a closed path using a random set of points
  int npoints = 5;
  path = hmap::Path(npoints, seed, bbox);
  path.set_closed(true);

  path.reorder_nns(); // reorder points to get a better look

  hmap::Array z2 = hmap::Array(shape);
  path.to_array(z2, bbox);

  // fractalize
  int   iterations = 4;
  float sigma = 0.3f;

  path.resample_uniform(); // to ensure a "uniform" output
  path = hmap::fractalize(path, iterations, seed, sigma);

  hmap::Array z3 = hmap::Array(shape);
  path.to_array(z3, bbox);

  hmap::export_banner_png("ex_path.png", {z1, z2, z3}, hmap::Cmap::INFERNO);
}
