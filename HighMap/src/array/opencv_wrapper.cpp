/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/core/utility.hpp>

#include "highmap/adapters/opencv.hpp"
#include "highmap/array.hpp"
#include "highmap/range.hpp"

namespace hmap
{

cv::Mat to_cv_mat(Array &array)
{
  return cv::Mat(array.shape.y, array.shape.x, CV_32FC1, array.vector.data());
}

// Helper function to convert OpenCV matrix to Array
template <typename T>
void convert_mat_to_array(const cv::Mat &mat, Array &array, bool flip_j)
{
  for (int j = 0; j < array.shape.y; ++j)
  {
    int jj = flip_j ? (array.shape.y - 1 - j) : j;

    for (int i = 0; i < array.shape.x; ++i)
      array(i, jj) = static_cast<float>(mat.at<T>(j, i));
  }
}

Array cv_mat_to_array(const cv::Mat &mat, bool remap_values, bool flip_j)
{
  glm::ivec2 shape = {mat.cols, mat.rows};
  Array      array(shape);

  switch (mat.type())
  {
  case CV_8U: convert_mat_to_array<uint8_t>(mat, array, flip_j); break;
  case CV_8S: convert_mat_to_array<int8_t>(mat, array, flip_j); break;
  case CV_16U: convert_mat_to_array<uint16_t>(mat, array, flip_j); break;
  case CV_16S: convert_mat_to_array<int16_t>(mat, array, flip_j); break;
  case CV_32S: convert_mat_to_array<int>(mat, array, flip_j); break;
  case CV_32F: convert_mat_to_array<float>(mat, array, flip_j); break;
  case CV_64F: convert_mat_to_array<double>(mat, array, flip_j); break;
  default: throw std::invalid_argument("Unsupported matrix type.");
  }

  if (remap_values) remap(array);

  return array;
}

std::string get_opencv_build_information()
{
  return cv::getBuildInformation();
}

} // namespace hmap
