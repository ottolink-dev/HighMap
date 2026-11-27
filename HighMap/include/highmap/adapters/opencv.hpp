/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file opencv_bridge.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief OpenCV cv::Mat <-> hmap::Array conversion helpers.
 * @copyright Copyright (c) 2023
 */
#pragma once

#include <string>

#include <opencv2/core.hpp>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Converts an OpenCV `cv::Mat` to a 2D `Array` with optional value
 * scaling to [0, 1].
 */
Array cv_mat_to_array(const cv::Mat &mat,
                      bool           remap_values = true,
                      bool           flip_j = false);

/**
 * @brief Returns OpenCV build information string.
 */
std::string get_opencv_build_information();

/**
 * @brief Converts a 2D `Array` to a non-owning OpenCV `cv::Mat` (CV_32FC1).
 */
cv::Mat to_cv_mat(Array &array);

} // namespace hmap
