#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 256};
  glm::vec2  res = {4.f, 2.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z);

  z.infos();

  cv::Mat mat = hmap::to_cv_mat(z);

  // accepts various matrix CV types and if multiple channels are
  // available, only the first one is used

  // mat.convertTo(mat, CV_8UC1, 255);
  // mat.convertTo(mat, CV_16UC3, 65535);
  // mat.convertTo(mat, CV_64F);

  hmap::Array zm = hmap::cv_mat_to_array(mat);

  hmap::export_banner_png("ex_from_cv_mat.png", {z, zm}, hmap::Cmap::MAGMA);
}
