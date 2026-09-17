#include <iostream>

#include "highmap.hpp"

int main(void)
{
  // generate points
  int        npoints = 7;
  int        seed = 1;
  glm::vec4  bbox = {-1.f, 0.f, 0.5f, 1.5f};
  hmap::Path path = hmap::Path(npoints, seed, bbox);
  path.reorder_nns();

  // interpolate
  std::vector<hmap::Point> pts(path.begin(), path.end());

  hmap::InterpolatorCurve fb = hmap::InterpolatorCurve(
      pts,
      hmap::InterpolationMethodCurve::BEZIER);

  hmap::InterpolatorCurve fs = hmap::InterpolatorCurve(
      pts,
      hmap::InterpolationMethodCurve::BSPLINE);

  hmap::InterpolatorCurve fc = hmap::InterpolatorCurve(
      pts,
      hmap::InterpolationMethodCurve::CATMULLROM);

  hmap::InterpolatorCurve fd = hmap::InterpolatorCurve(
      pts,
      hmap::InterpolationMethodCurve::DECASTELJAU);

  int                npts = 200;
  std::vector<float> t = hmap::linspace(0.f, 1.f, npts);

  // plots
  hmap::Array z = hmap::Array(glm::ivec2(512, 512));

  path.to_array(z, bbox);
  hmap::Path(fb(t)).to_array(z, bbox);
  hmap::Path(fs(t)).to_array(z, bbox);
  hmap::Path(fc(t)).to_array(z, bbox);
  hmap::Path(fd(t)).to_array(z, bbox);

  z.to_png("ex_interpolate_curve.png", hmap::Cmap::INFERNO);

  hmap::Path(fd(t)).to_png("out.png");
}
