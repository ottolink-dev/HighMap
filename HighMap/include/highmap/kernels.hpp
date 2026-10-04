/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file kernels.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for kernel functions and utilities.
 *
 * This header file declares functions for generating various types of kernel
 * functions used in statistical and computational applications. These kernels
 * include common types such as biweight, cubic pulse, Lorentzian, and smooth
 * cosine, among others.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <functional>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Enumeration for different kernel functions used in various algorithms.
 *
 * This enumeration defines the types of kernel functions that can be used for
 * smoothing, interpolation, and other operations requiring a kernel. Each
 * kernel type represents a specific mathematical function used to weight data
 * points based on their distance from a central point.
 *
 * These kernels are used in algorithms that require weighting functions, such
 * as kernel density estimation, interpolation, and data smoothing.
 */
enum KernelType : int
{
  BIWEIGHT,      ///< Biweight kernel function
  CUBIC_PULSE,   ///< Cubic pulse kernel function
  CONE,          ///< Cone kernel function
  CONE_SMOOTH,   ///< Smooth cone kernel function
  DISK,          ///< Disk kernel function
  LORENTZIAN,    ///< Lorentzian (Cauchy) kernel function
  SMOOTH_COSINE, ///< Smooth cosine kernel function
  SQUARE,        ///< Square kernel function
  TRICUBE,       ///< Tricube kernel function
  CUPOLA         ///< Cupola (reverse) kernel function
};

/**
 * @brief Generate a kernel of the specified type.
 *
 * This function creates a kernel based on the given shape and kernel type. It
 * supports various types of kernels including biweight, cubic pulse, cone, and
 * more. The resulting kernel is returned as an array.
 *
 * @param  shape       The dimensions of the kernel to be generated.
 * @param  kernel_type The type of kernel to generate (e.g., BIWEIGHT,
 *                     CUBIC_PULSE, etc.).
 * @return             Array The generated kernel array.
 */
Array get_kernel(glm::ivec2 shape, KernelType kernel_type);

// --- KERNELS

/**
 * @brief Generates a biweight kernel array.
 *
 * This function creates a biweight kernel, which is a type of kernel function
 * used in statistics for smoothing and interpolation.
 *
 * For more details on the biweight kernel, refer to the Wikipedia page:
 * https://en.wikipedia.org/wiki/Kernel_%28statistics%29.
 *
 * @param  shape Array shape specifying the dimensions of the kernel.
 * @return       Array A new array containing the biweight kernel.
 */
Array biweight(glm::ivec2 shape);

/**
 * @brief Generates a Blackman window array with the specified shape.
 *
 * The Blackman window is commonly used in signal processing for smoothing or
 * tapering. This function creates a 2D array of values representing the
 * Blackman window, calculated based on the given shape dimensions.
 *
 * @param  shape A 2D vector representing the dimensions of the array.
 *
 * @return       An `Array` object containing the Blackman window values of the
 *               specified shape.
 */
Array blackman(glm::ivec2 shape);

/**
 * @brief Generates a cone-shaped kernel array.
 *
 * This function creates a cone-shaped kernel. The maximum value of the kernel
 * is 1, and the shape of the kernel is determined by the provided array
 * dimensions.
 *
 * @param  shape Array shape specifying the dimensions of the kernel.
 * @return       Array A new array containing the cone-shaped kernel.
 */
Array cone(glm::ivec2 shape);

/**
 * @brief Generates a cone-shaped kernel with specified height and talus.
 *
 * This function creates a cone-shaped kernel with a given height and talus. The
 * output array's shape is adjusted accordingly based on the specified height
 * and talus.
 *
 * @param  height Height of the cone. This determines the peak value of the
 *                kernel.
 * @param  talus  The slope of the cone, which affects the rate at which the
 *                kernel value decreases.
 * @return        Array A new array containing the cone-shaped kernel with the
 *                specified height and talus.
 */
Array cone_talus(float height, float talus);

/**
 * @brief Generates a cone-shaped kernel with a smooth landing.
 *
 * This function creates a cone-shaped kernel where the kernel smoothly
 * transitions to zero at the bottom, resulting in a smooth derivative at the
 * cone base.
 *
 * @param  shape Array shape specifying the dimensions of the kernel.
 * @return       Array A new array containing the smooth cone-shaped kernel.
 */
Array cone_smooth(glm::ivec2 shape);

/**
 * @brief Generates a cubic pulse kernel array.
 *
 * This function creates a cubic pulse kernel, which is a type of kernel
 * function characterized by its cubic shape.
 *
 * @param  shape Array shape specifying the dimensions of the kernel.
 * @return       Array A new array containing the cubic pulse kernel.
 */
Array cubic_pulse(glm::ivec2 shape);

/**
 * @brief Generates a 1D cubic pulse kernel.
 *
 * @param  nk The number of samples in the 1D kernel.
 * @return    std::vector<float> A vector containing the cubic pulse kernel
 *            values.
 */
std::vector<float> cubic_pulse_1d(int nk);

/**
 * @brief Generates a "directional" cubic pulse kernel.
 *
 * This function creates a cubic pulse kernel with directional properties. The
 * kernel's shape can be adjusted based on the provided angle, aspect ratio, and
 * anisotropy. This allows for creating kernels that are elongated or compressed
 * in specific directions, making them suitable for applications requiring
 * directional smoothing or filtering.
 *
 * @param  shape        Array shape specifying the dimensions of the kernel.
 * @param  angle        Angle (in degrees) defining the direction of the pulse
 *                      elongation.
 * @param  aspect_ratio Pulse aspect ratio, which controls the kernel's
 *                      elongation or compression.
 * @param  anisotropy   Pulse width ratio between upstream and downstream sides
 *                      of the pulse.
 * @return              Array A new array containing the directional cubic pulse
 *                      kernel.
 */
Array cubic_pulse_directional(glm::ivec2 shape,
                              float      angle,
                              float      aspect_ratio,
                              float      anisotropy);

/**
 * @brief Generates a truncated cubic pulse kernel.
 *
 * This function creates a cubic pulse kernel that is truncated based on the
 * specified slant ratio and angle. The truncation slope determines how the
 * kernel's values decrease as they approach the edge of the kernel, creating a
 * "cut-off" effect.
 *
 * @param  shape       Array shape specifying the dimensions of the kernel.
 * @param  slant_ratio Truncation slope that affects the kernel's cut-off
 *                     behavior.
 * @param  angle       Angle (in degrees) that defines the orientation of the
 *                     truncation.
 * @return             Array A new array containing the truncated cubic pulse
 *                     kernel.
 *
 * **Example**
 * @include ex_cubic_pulse_truncated.cpp
 *
 * **Result**
 * @image html ex_cubic_pulse_truncated.png
 */
Array cubic_pulse_truncated(glm::ivec2 shape, float slant_ratio, float angle);

/**
 * @brief Generate a 2D radial cupola kernel.
 *
 * Values are computed from the normalized radial distance to the center:
 * constant 1 inside radius @p rc, then smoothly decreasing to 0 at radius 1.
 *
 * @param  shape Size of the output array (width, height).
 * @param  rc    Radius of the inner flat region (normalized, < 1).
 * @return       Generated cupola-shaped Array.
 */
Array cupola(glm::ivec2 shape, float rc = 0.5f);

/**
 * @brief Generates a disk-shaped kernel footprint.
 *
 * This function creates a disk-shaped kernel, which is a circular kernel where
 * the values are radially symmetric around the center.
 *
 * @param  shape Array shape specifying the dimensions of the kernel.
 * @return       Array A new array containing the disk-shaped kernel.
 */
Array disk(glm::ivec2 shape);

/**
 * @brief Generates a smooth, disk-shaped kernel footprint with soft edges.
 *
 * This function creates a 2D, disk-shaped kernel with radially symmetric values
 * centered in the array. Values inside the cutoff radius `r_cutoff` are set to
 * 1.0. Between `r_cutoff` and the normalized radius 1.0, the values smoothly
 * fall off to zero using a parabolic profile to ensure soft transitions at the
 * disk boundary.
 *
 * The generated array is useful for constructing smooth low-pass filters or
 * spatial masks with circular support.
 *
 * @param  shape    The dimensions of the output array (width, height).
 * @param  r_cutoff Normalized cutoff radius in the range [0.0, 1.0]. Inside
 *                  this radius, the kernel has a flat value of 1.0; outside, it
 *                  decays smoothly to 0.0.
 *
 * @return          A 2D array of size `shape` containing the smooth disk-shaped
 *                  kernel.
 *
 * @note Coordinates are normalized such that the radius 1.0 corresponds to the
 * maximum possible distance from the center in either dimension.
 *
 * @see             disk()
 */
Array disk_smooth(glm::ivec2 shape, float r_cutoff = 0.9f);

/**
 * @brief Generates a Gabor kernel of the specified shape.
 *
 * This function creates a Gabor kernel, which is a type of sinusoidal kernel
 * modulated by a Gaussian envelope. The Gabor kernel is useful for texture
 * analysis and feature extraction. The kernel's frequency and orientation are
 * defined by the wavenumber and angle parameters, respectively.
 *
 * @param  shape Array shape specifying the dimensions of the kernel.
 * @param  kw    Kernel wavenumber, which determines the frequency of the
 *               sinusoidal component.
 * @param  angle Kernel angle (in degrees) that defines the orientation of the
 *               sinusoidal component.
 * @return       Array A new array containing the Gabor kernel.
 *
 * **Example**
 * @include ex_gabor.cpp
 *
 * **Result**
 * @image html ex_gabor.png
 */
Array gabor(glm::ivec2 shape,
            float      kw,
            float      angle,
            bool       quad_phase_shift = false);

/**
 * @brief Generates a modified dune-like Gabor kernel.
 *
 * This function creates a modified Gabor kernel that has a dune-like profile.
 * The dune profile modifies the standard Gabor kernel to include a top and foot
 * for the dune, giving it a specific shape. The `xtop` and `xbottom` parameters
 * control the relative locations of the top and foot of the dune profile,
 * respectively.
 *
 * @param  shape   Array shape specifying the dimensions of the kernel.
 * @param  kw      Kernel wavenumber, which determines the frequency of the
 *                 sinusoidal component.
 * @param  angle   Kernel angle (in degrees) that defines the orientation of the
 *                 sinusoidal component.
 * @param  xtop    Relative location of the top of the dune profile (in [0, 1]).
 * @param  xbottom Relative location of the foot of the dune profile (in [0,
 *                 1]).
 * @return         Array A new array containing the modified dune-like Gabor
 *                 kernel.
 *
 * **Example**
 * @include ex_gabor_dune.cpp
 *
 * **Result**
 * @image html ex_gabor_dune.png
 */
Array gabor_dune(glm::ivec2 shape,
                 float      kw,
                 float      angle,
                 float      xtop,
                 float      xbottom);

/**
 * @brief Generates a Hann window array with the specified shape.
 *
 * The Hann window is commonly used in signal processing for smoothing or
 * tapering. This function creates a 2D array of values representing the Hann
 * window, calculated based on the given shape dimensions.
 *
 * @param  shape A 2D vector representing the dimensions of the array.
 *
 * @return       An `Array` object containing the Hann window values of the
 *               specified shape.
 */
Array hann(glm::ivec2 shape);

/**
 * @brief Generate a Lorentzian kernel.
 *
 * The Lorentzian kernel is characterized by its peak and tails that decay with
 * the distance from the center. The width of the kernel is determined using the
 * provided footprint threshold.
 *
 * @param  shape               The dimensions of the kernel.
 * @param  footprint_threshold Determines the width of the Lorentzian kernel.
 *                             Default is 0.1.
 * @return                     Array The Lorentzian kernel array.
 */
Array lorentzian(glm::ivec2 shape, float footprint_threshold = 0.1f);

/**
 * @brief Generate a modified Lorentzian kernel with compact support.
 *
 * This version of the Lorentzian kernel is modified to have compact support,
 * meaning it is non-zero only within a certain range. This modification is
 * useful for specific applications requiring compact kernels.
 *
 * @param  shape The dimensions of the kernel.
 * @return       Array The modified Lorentzian kernel with compact support.
 */
Array lorentzian_compact(glm::ivec2 shape);

/**
 * @brief Generates a radial sinc function array with the specified shape and
 * wave number.
 *
 * This function creates a 2D array of values representing the radial sinc
 * function, calculated based on the given shape and wave number. The sinc
 * function is defined as `sin(x) / x` and is useful in signal processing and
 * analysis.
 *
 * @param  shape A 2D vector representing the dimensions of the array.
 * @param  kw    The wave number used to scale the radial sinc function.
 *
 * @return       An `Array` object containing the radial sinc function values of
 *               the specified shape.
 */
Array sinc_radial(glm::ivec2 shape, float kw);

/**
 * @brief Generates a separable sinc function array with the specified shape and
 * wave number.
 *
 * This function creates a 2D array of values representing the separable sinc
 * function, calculated as the product of 1D sinc functions along each dimension
 * of the array. The sinc function is defined as `sin(x) / x` and is useful in
 * signal processing and analysis.
 *
 * @param  shape A 2D vector representing the dimensions of the array.
 * @param  kw    The wave number used to scale the separable sinc function.
 *
 * @return       An `Array` object containing the separable sinc function values
 *               of the specified shape.
 */
Array sinc_separable(glm::ivec2 shape, float kw);

/**
 * @brief Generate a smooth cosine kernel.
 *
 * The smooth cosine kernel is a type of kernel function that smoothly
 * transitions between values, based on a cosine function.
 *
 * @param  shape The dimensions of the kernel.
 * @return       Array The smooth cosine kernel array.
 */
Array smooth_cosine(glm::ivec2 shape);

/**
 * @brief Generate a square-shaped kernel.
 *
 * The square kernel is a simple, isotropic kernel where the values are constant
 * within the square region. The maximum value within the kernel is 1.
 *
 * @param  shape The dimensions of the kernel.
 * @return       Array The square-shaped kernel array.
 */
Array square(glm::ivec2 shape);

/**
 * @brief Generate a tricube kernel.
 *
 * The tricube kernel is a type of kernel function that is often used in
 * statistical applications. It has a specific shape defined by the tricube
 * polynomial.
 *
 * @param  shape The dimensions of the kernel.
 * @return       Array The tricube kernel array.
 *
 * **Reference**
 * See https://en.wikipedia.org/wiki/Kernel_%28statistics%29 for details on the
 * tricube kernel.
 */
Array tricube(glm::ivec2 shape);

} // namespace hmap
