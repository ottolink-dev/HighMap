/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file roads.hpp
 * @author  Otto Link (otto.link.bv@gmail.com)
 * @brief
 * @version 0.1
 * @date 2023-06-26
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <functional>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Define on which part of the pyramid component the transform function
 * should be applied.
 */
enum pyramid_transform_support : int
{
  FULL,         ///< Complete field (sum of low and high pass)
  LOWPASS_ONLY, ///< Low pass component only
  HIGHPASS_ONLY ///< High pass component only
};

/**
 * @brief Pyramid decomposition class, to handle low-pass pyramids (like
 * Laplacian pyramid).
 */
class PyramidDecomposition
{
public:
  /**
   * @brief Number of levels in the pyramid.
   */
  int nlevels;

  /**
   * @brief Residual field (low-pass component) a the coarsest level.
   */
  Array residual;

  /**
   * @brief High-pass component for each level
   */
  std::vector<Array> components = {};

  /**
   * @brief Reference to the low-pass filter function.
   */
  std::function<Array(const Array &)> low_pass_filter_function;

  /**
   * @brief Construct a new Pyramid Decomposition object.
   *
   * @param array   Reference to the input array.
   * @param nlevels Number of levels (if set to a null or negative value, the
   *                maximum number of levels is taken minus the number
   *                provided).
   *
   * **Example**
   * @include ex_pyramid_decomposition.cpp
   *
   * **Result**
   * @image html ex_pyramid_decomposition.png
   */
  PyramidDecomposition(Array &array, int nlevels);

  /**
   * @brief Generate the pyramid decomposition.
   */
  void decompose();

  /**
   * @brief Reconstruct a field (ideally identical to the input field) based on
   * the pyramid decomposition.
   * @return Array Reconstructed field.
   */
  Array reconstruct();

  /**
   * @brief Export pyramid as png image file.
   *
   * @param fname       File name.
   * @param cmap        Colormap (@see cmap).
   * @param hillshading Activate hillshading.
   */
  void to_png(std::string fname, int cmap, bool hillshading = false);

  /**
   * @brief Apply a transformation to each component of the pyramid and return
   * the field constructed based on the modified pyramid.
   *
   * @param  function      Reference to the function applied to the pyramid
   *                       components.
   * @param  support       Function support, should it be applied to the lowpass
   *                       components only, the highpass only the full field.
   * @param  level_weights Weight in [0, 1] for each level (the resulting
   *                       component is lerp between no transform and transform
   *                       according to this weight).
   * @return               Array Resulting array.
   *
   * **Example**
   * @include ex_pyramid_transform.cpp
   *
   * **Result**
   * @image html ex_pyramid_transform.png
   */
  Array transform(
      std::function<Array(const Array &, const int current_level)> function,
      int                                                          support = 0,
      std::vector<float> level_weights = {},
      int                finest_level = 0);

private:
  /**
   * @brief Reference to the input array.
   */
  Array *p_array = nullptr;
};

} // namespace hmap