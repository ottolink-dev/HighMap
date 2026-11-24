/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>

#include "highmap/adapters/opencv.hpp"
#include "highmap/array.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

Array read_to_array(const std::string &fname, bool flip_j, bool remap)
{
  cv::Mat mat = cv::imread(fname, cv::IMREAD_GRAYSCALE | cv::IMREAD_ANYDEPTH);

  if (mat.data == nullptr)
  {
    hmap::log::error("error while reading the image file: {}", fname);
    return Array();
  }
  else
  {
    return cv_mat_to_array(mat, remap, flip_j);
  }
}

} // namespace hmap
