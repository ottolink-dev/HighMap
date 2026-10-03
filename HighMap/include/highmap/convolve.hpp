/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file convolve.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Provides functions for convolution operations on arrays.
 *
 * This header file defines functions for performing convolution operations on
 * 1D and 2D arrays using various kernels. The operations include 1D
 * convolutions along different directions, 2D convolutions with symmetry
 * boundary conditions, truncated 2D convolutions, and approximate convolutions
 * using Singular Value Decomposition (SVD).
 *
 * @copyright Copyright (c) 2023 Otto Link
 */

#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Return the convolution product of the array with a 1D kernel along the
 * 'i' direction.
 *
 * This function applies a 1D convolution to the input array using the provided
 * kernel, considering the convolution in the row direction ('i' direction).
 *
 * @param  array  Input array to be convolved.
 * @param  kernel 1D kernel to be used for the convolution.
 * @return        Array Resulting array after applying the 1D convolution.
 *
 * **Example**
 * @include ex_convolve1d_ij.cpp
 *
 * **Result**
 * @image html ex_convolve1d_ij.png
 *
 * @see           {@link convolve1d_j}
 */
Array convolve1d_i(const Array &array, const std::vector<float> &kernel);

/**
 * @brief Return the convolution product of the array with a 1D kernel along the
 * 'j' direction.
 *
 * This function applies a 1D convolution to the input array using the provided
 * kernel, considering the convolution in the column direction ('j' direction).
 *
 * @param  array  Input array to be convolved.
 * @param  kernel 1D kernel to be used for the convolution.
 * @return        Array Resulting array after applying the 1D convolution.
 *
 * **Example**
 * @include ex_convolve1d_ij.cpp
 *
 * **Result**
 * @image html ex_convolve1d_ij.png
 *
 * @see           {@link convolve1d_i}
 */
Array convolve1d_j(const Array &array, const std::vector<float> &kernel);

/**
 * @brief Return the convolution product of the array with a given 2D kernel.
 *
 * This function performs a 2D convolution on the input array using the
 * specified 2D kernel. The output array has the same shape as the input array,
 * and symmetry boundary conditions are used.
 *
 * @param  array  Input array to be convolved.
 * @param  kernel 2D kernel to be used for the convolution.
 * @return        Array Resulting array after applying the 2D convolution.
 *
 * **Example**
 * @include ex_convolve2d_svd.cpp
 */
Array convolve2d(const Array &array, const Array &kernel);

/**
 * @brief Return the convolution product of the array with a given 2D kernel,
 * with a truncated output size.
 *
 * This function performs a 2D convolution on the input array using the
 * specified 2D kernel. The resulting output array is smaller than the input
 * array by the dimensions of the kernel.
 *
 * @param  array  Input array to be convolved.
 * @param  kernel 2D kernel to be used for the convolution.
 * @return        Array Resulting array after applying the truncated 2D
 *                convolution. Shape: {array.shape[0] - kernel.shape[0],
 *                array.shape[1] -
 *                kernel.shape[1]}.
 */
Array convolve2d_truncated(const Array &array, const Array &kernel);

/**
 * @brief Return the approximate convolution product of the array with a
 * Singular Value Decomposition (SVD) of a kernel.
 *
 * This function approximates the convolution of the input array with a 2D
 * kernel by using its Singular Value Decomposition (SVD). The approximation is
 * based on the first `rank` singular values and vectors.
 *
 * @param  z      Input array to be convolved.
 * @param  kernel 2D kernel to be used for the convolution.
 * @param  rank   Approximation rank: number of singular values/vectors used for
 *                the approximation.
 * @return        Array Resulting array after applying the SVD-based convolution
 *                approximation.
 *
 * **Example**
 * @include ex_convolve2d_svd.cpp
 *
 * **Result**
 * @image html ex_convolve2d_svd.png
 */
Array convolve2d_svd(const Array &z, const Array &kernel, int rank = 3);

/**
 * @brief Return the approximate convolution product of the array with a
 * Singular Value Decomposition (SVD) of a kernel combined with kernel
 * rotations.
 *
 * This function approximates the convolution of the input array with a 2D
 * kernel using its SVD, combined with multiple rotations of the kernel. This
 * can enhance the convolution approximation by considering different
 * orientations of the kernel.
 *
 * @param  z           Input array to be convolved.
 * @param  kernel      2D kernel to be used for the convolution.
 * @param  rank        Approximation rank: number of singular values/vectors
 *                     used for the approximation.
 * @param  n_rotations Number of kernel rotations to be considered.
 * @param  seed        Random seed number for kernel rotations.
 * @return             Array Resulting array after applying the SVD-based
 *                     convolution approximation with rotations.
 *
 * **Example**
 * @include ex_convolve2d_svd.cpp
 *
 * **Result**
 * @image html ex_convolve2d_svd.png
 */
Array convolve2d_svd_rotated_kernel(const Array  &z,
                                    const Array  &kernel,
                                    int           rank = 3,
                                    int           n_rotations = 6,
                                    std::uint32_t seed = 1);

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Apply a sparse max-convolution.
 *
 * This operator propagates the non-zero values of the input array through a
 * spatial kernel using a max-reduction instead of a summation. For each
 * non-zero value of the input array, the kernel is centered at the sample
 * position and contributes to the output through a weighted maximum.
 *
 * The operation is conceptually similar to a sparse morphological dilation or a
 * max-convolution.
 *
 * @param  array                   Input sparse scalar field.
 * @param  kernel                  Propagation kernel.
 * @param  scale_kernel_with_value If true, the effective kernel size is scaled
 *                                 according to the local input value. Larger
 *                                 values propagate over a larger neighborhood,
 *                                 while smaller values use a reduced kernel
 *                                 support. If false, the full kernel size is
 *                                 used everywhere.
 * @param  k_smoothmax             Smooth maximum factor.
 *
 * @return                         Sparse max-convolved array.
 *
 * **Example**
 * @include ex_sparse_max_convolution.cpp
 *
 * **Result**
 * @image html ex_sparse_max_convolution.png
 */
Array sparse_max_convolution(const Array &array,
                             const Array &kernel,
                             bool         scale_kernel_with_value = true,
                             float        k_smoothmax = 0.f);

Array sparse_max_convolution(const Array &array_values,
                             const Array &array_sizes,
                             const Array &kernel,
                             float        k_smoothmax = 0.f);

} // namespace hmap::gpu