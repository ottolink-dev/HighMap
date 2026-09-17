#include "highmap.hpp"

int main(void)
{
  std::vector<float> x = {0.f, 0.1f, 0.5f, 0.7f};
  std::vector<float> y = {0.2f, 0.1f, 0.5f, 0.6f};

  // generate 4 points
  glm::vec4  bbox = {-1.f, 0.f, 0.5f, 1.5f};
  hmap::Path path = hmap::Path(4, 2, bbox);
  path.reorder_nns();

  path.print();

  hmap::Point p0 = path[0];
  hmap::Point p1 = path[1];
  hmap::Point p2 = path[2];
  hmap::Point p3 = path[3];

  // interpolate
  int                npts = 50;
  std::vector<float> t = hmap::linspace(0.f, 1.f, npts);

  std::vector<hmap::Point> points_bezier = {};
  std::vector<hmap::Point> points_bspline = {};
  std::vector<hmap::Point> points_catmullrom = {};

  for (auto &t_ : t)
  {
    points_bezier.push_back(hmap::interp_bezier(p0, p1, p2, p3, t_));
    points_bspline.push_back(hmap::interp_bspline(p0, p1, p2, p3, t_));
    points_catmullrom.push_back(hmap::interp_catmullrom(p0, p1, p2, p3, t_));
  }

  // plots
  hmap::Array z = hmap::Array(glm::ivec2(512, 512));
  path.to_array(z, bbox);
  hmap::Path(points_bezier).to_array(z, bbox);
  hmap::Path(points_bspline).to_array(z, bbox);
  hmap::Path(points_catmullrom).to_array(z, bbox);
  z.to_png("ex_point_interp.png", hmap::Cmap::INFERNO);
}
