#include <opencv2/opencv.hpp>

#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 256};
  glm::vec2  res = {4.f, 2.f};
  int        seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);
  hmap::remap(z, 0.f, 1.f);

  z.to_png("out.png", hmap::Cmap::GRAY);

  z.infos();

  cv::Mat img = hmap::to_cv_mat(z);
  cv::cvtColor(img, img, cv::COLOR_GRAY2BGR);

  std::vector<int> codec_params = {cv::IMWRITE_EXR_TYPE,
                                   cv::IMWRITE_EXR_TYPE_FLOAT,
                                   cv::IMWRITE_EXR_COMPRESSION,
                                   cv::IMWRITE_EXR_COMPRESSION_NO};

  cv::imwrite("out.exr", img, codec_params);

  // set compression to cv::IMWRITE_TIFF_COMPRESSION_LZW (apparently
  // not available in openCV public header?)
  codec_params = {cv::IMWRITE_TIFF_COMPRESSION, 5};

  cv::imwrite("out.tiff", img, codec_params);

  img.convertTo(img, CV_16U, 65535);
  cv::imwrite("out_cv.png", img);

  cv::imshow("example", img);
  cv::waitKey(0);
}
