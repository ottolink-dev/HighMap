#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  int        seed = 0;

  hmap::BBox2D bbox = {1.f, 2.f, -0.5f, 0.5f};
  hmap::Path   path = hmap::Path(10,
                               seed,
                               hmap::BBox2D({1.2f, 1.8f, -0.3f, 0.3f}));
  path.reorder_nns();

  const int npoints = 100;

  const std::vector<hmap::InterpolationMethod1D> methods = {
      hmap::InterpolationMethod1D::AKIMA,
      hmap::InterpolationMethod1D::AKIMA_PERIODIC,
      hmap::InterpolationMethod1D::CUBIC,
      hmap::InterpolationMethod1D::CUBIC_PERIODIC,
      hmap::InterpolationMethod1D::LINEAR,
      // hmap::InterpolationMethod1D::POLYNOMIAL, // severe overshoot
      // hmap::InterpolationMethod1D::STEFFEN // cannot be used
  };

  std::vector<hmap::Array> arrays;

  // --- Open path

  path.set_closed(false);

  // original path
  arrays.clear();
  arrays.push_back(path.to_array(shape, bbox));

  for (const auto method : methods)
  {
    hmap::Path path_copy = path;
    path_copy.resample_interp(npoints, method);

    // plot to an array
    arrays.push_back(path_copy.to_array(shape, bbox));
  }

  hmap::export_banner_png("ex_path_resample0.png", arrays, hmap::Cmap::INFERNO);

  // --- Closed path

  path.set_closed(true);

  // original path
  arrays.clear();
  arrays.push_back(path.to_array(shape, bbox));

  for (const auto method : methods)
  {
    hmap::Path path_copy = path;
    path_copy.resample_interp(npoints, method);

    // plot to an array
    arrays.push_back(path_copy.to_array(shape, bbox));
  }

  hmap::export_banner_png("ex_path_resample1.png", arrays, hmap::Cmap::INFERNO);

  // --- Wrappers

  path.set_closed(false);

  arrays.clear();
  arrays.push_back(path.to_array(shape, bbox));

  {
    hmap::Path path_copy = path;
    path_copy.resample_uniform();
    arrays.push_back(path_copy.to_array(shape, bbox));
  }

  {
    hmap::Path path_copy = path;
    path_copy.resample_by_spacing(0.01f);
    arrays.push_back(path_copy.to_array(shape, bbox));
  }

  hmap::export_banner_png("ex_path_resample2.png", arrays, hmap::Cmap::INFERNO);
}
