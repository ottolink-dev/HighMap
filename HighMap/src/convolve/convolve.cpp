/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cmath>
#include <vector>

#include <opencv2/core.hpp>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/convolve.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

Array convolve1d_i(const Array &array, const std::vector<float> &kernel)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_non_empty(kernel, "Kernel")) return Array(array.shape);

  Array     array_out = Array(array.shape);
  const int nk = (int)kernel.size();
  const int i1 = nk / 2;

  for (int p = 0; p < nk; p++)
  {
    for (int i = 0; i < array.shape.x; i++)
    {
      const int ii = std::clamp(i + p - i1, 0, array.shape.x - 1);
      for (int j = 0; j < array.shape.y; j++)
        array_out(i, j) += array(ii, j) * kernel[p];
    }
  }
  return array_out;
}

Array convolve1d_j(const Array &array, const std::vector<float> &kernel)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_non_empty(kernel, "Kernel")) return Array(array.shape);

  Array     array_out = Array(array.shape);
  const int nk = (int)kernel.size();
  const int j1 = nk / 2;

  for (int p = 0; p < nk; p++)
  {
    for (int i = 0; i < array.shape.x; i++)
    {
      for (int j = 0; j < array.shape.y; j++)
      {
        const int jj = std::clamp(j + p - j1, 0, array.shape.y - 1);
        array_out(i, j) += array(i, jj) * kernel[p];
      }
    }
  }
  return array_out;
}

Array convolve2d(const Array &array, const Array &kernel)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_non_empty(kernel)) return Array(array.shape);

  const int i1 = (int)ceil(0.5f * (float)kernel.shape.x);
  const int i2 = kernel.shape.x - i1;
  const int j1 = (int)ceil(0.5f * (float)kernel.shape.y);
  const int j2 = kernel.shape.y - j1;

  Array array_buffered = generate_buffered_array(array, {i1, i2, j1, j2});
  Array array_out = convolve2d_truncated(array_buffered, kernel);

  return array_out;
}

Array convolve2d_truncated(const Array &array, const Array &kernel)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_non_empty(kernel)) return Array();

  if (array.shape.x <= kernel.shape.x || array.shape.y <= kernel.shape.y)
  {
    hmap::log::warn(
        "Array shape ({}, {}) must be strictly larger than kernel shape ({}, "
        "{}) for truncated convolution",
        array.shape.x,
        array.shape.y,
        kernel.shape.x,
        kernel.shape.y);
    return Array();
  }

  const int out_x = array.shape.x - kernel.shape.x;
  const int out_y = array.shape.y - kernel.shape.y;
  Array     array_out = Array(glm::ivec2(out_x, out_y));

  // use direct spatial convolution for small kernels, FFT-based convolution
  // otherwise
  if (kernel.shape.x * kernel.shape.y < 64)
  {
    for (int j = 0; j < array_out.shape.y; j++)
      for (int i = 0; i < array_out.shape.x; i++)
        for (int q = 0; q < kernel.shape.y; q++)
          for (int p = 0; p < kernel.shape.x; p++)
            array_out(i, j) += array(i + p, j + q) * kernel(p, q);
  }
  else
  {
    // --- FFT-based 2D convolution using OpenCV DFT

    const int dft_rows = cv::getOptimalDFTSize(array.shape.y);
    const int dft_cols = cv::getOptimalDFTSize(array.shape.x);

    cv::Mat mat_a(array.shape.y,
                  array.shape.x,
                  CV_32F,
                  const_cast<float *>(array.vector.data()));
    cv::Mat padded_a;
    cv::copyMakeBorder(mat_a,
                       padded_a,
                       0,
                       dft_rows - array.shape.y,
                       0,
                       dft_cols - array.shape.x,
                       cv::BORDER_CONSTANT,
                       cv::Scalar::all(0));

    cv::Mat mat_k(kernel.shape.y,
                  kernel.shape.x,
                  CV_32F,
                  const_cast<float *>(kernel.vector.data()));
    cv::Mat flipped_k;
    cv::flip(mat_k, flipped_k, -1);
    cv::Mat padded_k;
    cv::copyMakeBorder(flipped_k,
                       padded_k,
                       0,
                       dft_rows - kernel.shape.y,
                       0,
                       dft_cols - kernel.shape.x,
                       cv::BORDER_CONSTANT,
                       cv::Scalar::all(0));

    cv::Mat dft_a, dft_k, dft_out, conv_full;
    cv::dft(padded_a, dft_a, cv::DFT_COMPLEX_OUTPUT);
    cv::dft(padded_k, dft_k, cv::DFT_COMPLEX_OUTPUT);
    cv::mulSpectrums(dft_a, dft_k, dft_out, 0, false);
    cv::idft(dft_out, conv_full, cv::DFT_SCALE | cv::DFT_REAL_OUTPUT);

    cv::Rect roi(kernel.shape.x - 1, kernel.shape.y - 1, out_x, out_y);
    cv::Mat  valid_region = conv_full(roi);

    for (int j = 0; j < out_y; ++j)
    {
      const float *row_ptr = valid_region.ptr<float>(j);
      for (int i = 0; i < out_x; ++i)
        array_out(i, j) = row_ptr[i];
    }
  }

  return array_out;
}

} // namespace hmap
