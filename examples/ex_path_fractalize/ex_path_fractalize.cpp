#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 2;

  hmap::BBox2D bbox = {-1.f, 2.f, 0.f, 5.f};

  // --- generate a path using a random set of points

  int        npoints = 8;
  hmap::Path path = hmap::Path(npoints, seed, bbox);
  path.reorder_nns(); // reorder points to get a better look

  hmap::Array z1 = path.to_array(shape, bbox);

  // --- control function (supposed to be in [0, 1])

  hmap::Array z_control = hmap::slope(shape, 0.f, -1.f);
  hmap::remap(z_control);

  int   iterations = 6;
  float sigma = 0.3f;
  path.resample_uniform(); // to ensure a "uniform" output

  // --- fractalize, standard, bounded, and with control function

  int   orientation = 0;
  float persistence = 1.f;

  // 1. Standard unbounded fractalize
  hmap::Path pn = path;
  pn = hmap::fractalize(pn, iterations, seed, sigma);
  hmap::Array z3 = hmap::Array(shape);
  pn.to_array(z3, bbox);

  // 2. Bounded fractalize
  hmap::Path pb = path;
  pb = hmap::fractalize(pb,
                        iterations,
                        seed,
                        2.f * sigma, // push amplitude
                        orientation,
                        persistence,
                        nullptr,
                        bbox,
                        true);
  hmap::Array z5 = hmap::Array(shape);
  pb.to_array(z5, bbox);

  // 3. Control function modulated fractalize (closed path)
  hmap::Path pc = path;
  pc.set_closed(true);
  pc = hmap::fractalize(pc,
                        iterations,
                        seed,
                        sigma,
                        orientation,
                        persistence,
                        &z_control,
                        bbox);

  hmap::Array z4 = pc.to_array(shape, bbox);

  hmap::export_banner_png("ex_path_fractalize.png",
                          {z1, z3, z5, z_control, z4},
                          hmap::Cmap::INFERNO);
}
