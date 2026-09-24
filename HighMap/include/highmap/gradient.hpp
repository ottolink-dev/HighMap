/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file gradient.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Provides functions for calculating gradients and related operations on
 * arrays. This header file defines functions to compute various gradient
 * operations including gradient norms, gradient angles, and Laplacians for 2D
 * arrays. It supports different gradient filters such as Prewitt, Sobel, and
 * Scharr.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/texture.hpp"

namespace hmap
{

/**
 * @brief Compute the divergence of a 2D gradient field.
 *
 * Given two scalar fields representing the horizontal and vertical derivatives
 * of a height map:
 * \f[
 *     p(x, y) = \frac{\partial h}{\partial x}, \qquad q(x, y) = \frac{\partial
 * h}{\partial y}
 * \f] this function returns the divergence:
 * \f[
 *     f(x, y) = \frac{\partial p}{\partial x}(x, y)
 *             + \frac{\partial q}{\partial y}(x, y)
 * \f]
 *
 * The divergence is a key term for solving the Poisson equation
 * \f$\nabla^2 h = f\f$ when reconstructing a height field from gradients or
 * from a normal map.
 *
 * @param  dx Array of horizontal derivatives (\f$\partial h / \partial x\f$).
 *            Must have the same dimensions as @p dy.
 * @param  dy Array of vertical derivatives (\f$\partial h / \partial y\f$).
 *            Must have the same dimensions as @p dx.
 *
 * @return    A new Array containing the divergence \f$f(x, y)\f$ for each
 *            pixel. The returned array has the same size as the inputs.
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array divergence_from_gradients(const Array &dx, const Array &dy);

/**
 * @brief Compute the gradient of a 1D vector.
 *
 * This function calculates the gradient of a 1D vector by computing the
 * difference between successive elements.
 *
 * @param  v Input vector of floats.
 * @return   std::vector<float> Vector of gradient values.
 */
std::vector<float> gradient1d(const std::vector<float> &v);

/**
 * @brief Compute the polar angle of the gradient of a 2D array.
 *
 * This function calculates the gradient angle of a 2D array, which represents
 * the direction of the steepest ascent. The angle is measured in radians and
 * falls within the range [-π, π]. The direction can be inverted if the
 * `downward` flag is set to true.
 *
 * @param  array    Input 2D array.
 * @param  downward If true, invert the direction of the gradient.
 * @return          Array Gradient angles in radians.
 *
 * **Example**
 * @include ex_gradient_angle.cpp
 *
 * **Result**
 * @image html ex_gradient_angle.png
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_angle(const Array &array, bool downward = false);

/**
 * @brief Computes a smoothed gradient angle (aspect) field with circular
 * unwrapping.
 *
 * This function computes the gradient angle (aspect) of the input scalar field
 * and applies a circular smoothing to remove discontinuities caused by the
 * ±π wrap of the arctangent. The angle field is smoothed by converting the
 * angle to a unit vector representation on the circle, applying a local
 * smoothing kernel, and converting it back to an angle. This results in a
 * continuous aspect map that is ideal for visualization or further processing.
 *
 * @param  array    Input scalar field (e.g., a heightmap).
 * @param  ir       Radius of the local smoothing kernel applied to the angle
 *                  vectors.
 * @param  downward If set to `true`, angles are adjusted to point in the
 *                  downslope direction (useful for hydrology or flow
 *                  simulations). If `false`, angles point in the upslope
 *                  direction.
 *
 * @return          An Array containing the smoothed gradient angles in radians,
 *                  in the range [-π, π].
 *
 * @note
 * - Regions with very low gradient magnitude may produce unstable angle values;
 *   these regions should typically be masked or treated separately.
 * - This function is suited for preserving continuous angular fields when the
 * direction (not the orientation) of the gradient matters.
 *
 * **Example**
 * @include ex_gradient_angle.cpp
 *
 * **Result**
 * @image html ex_gradient_angle.png
 */
Array gradient_angle_circular_smoothing(const Array &array,
                                        int          ir,
                                        bool         downward = false);

/**
 * @brief Compute the gradient norm of a 2D array.
 *
 * This function calculates the gradient norm (magnitude) of a 2D array.
 * Optionally, the x and y gradients can be computed and stored in `p_dx` and
 * `p_dy` respectively.
 *
 * @param  array Input 2D array.
 * @param  p_dx  Pointer to an Array where the x-gradient will be stored
 *               (optional).
 * @param  p_dy  Pointer to an Array where the y-gradient will be stored
 *               (optional).
 * @return       Array Gradient norm (magnitude).
 *
 * **Example**
 * @include ex_gradient_norm.cpp
 *
 * **Result**
 * @image html ex_gradient_norm.png
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_norm(const Array &array,
                    Array       *p_dx = nullptr,
                    Array       *p_dy = nullptr);

/**
 * @brief Compute a smoothed gradient magnitude map.
 *
 * The input array is downsampled by a factor of @p ir, the gradient norm is
 * computed on the lower-resolution array, scaled back to the original units,
 * and finally upsampled to the original resolution using bicubic interpolation.
 *
 * @param  array Input array.
 * @param  ir    Downsampling factor (> 0).
 * @return       Gradient norm filtered at a larger spatial scale.
 *
 * **Example**
 * @include ex_gradient_norm_filtered.cpp
 *
 * **Result**
 * @image html ex_gradient_norm_filtered.png
 */
Array gradient_norm_filtered(const Array &array, int ir);

/**
 * @brief Compute the gradient norm of a 2D array using the Prewitt filter.
 *
 * This function calculates the gradient norm (magnitude) using the Prewitt
 * filter for gradient computation. Optionally, the x and y gradients can be
 * computed and stored in `p_dx` and `p_dy` respectively.
 *
 * @param  array Input 2D array.
 * @param  p_dx  Pointer to an Array where the x-gradient will be stored
 *               (optional).
 * @param  p_dy  Pointer to an Array where the y-gradient will be stored
 *               (optional).
 * @return       Array Gradient norm (magnitude).
 *
 * **Example**
 * @include ex_gradient_norm.cpp
 *
 * **Result**
 * @image html ex_gradient_norm.png
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_norm_prewitt(const Array &array,
                            Array       *p_dx = nullptr,
                            Array       *p_dy = nullptr);

/**
 * @brief Compute the gradient norm of a 2D array using the Scharr filter.
 *
 * This function calculates the gradient norm (magnitude) using the Scharr
 * filter for gradient computation. Optionally, the x and y gradients can be
 * computed and stored in `p_dx` and `p_dy` respectively.
 *
 * @param  array Input 2D array.
 * @param  p_dx  Pointer to an Array where the x-gradient will be stored
 *               (optional).
 * @param  p_dy  Pointer to an Array where the y-gradient will be stored
 *               (optional).
 * @return       Array Gradient norm (magnitude).
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_norm_scharr(const Array &array,
                           Array       *p_dx = nullptr,
                           Array       *p_dy = nullptr);

/**
 * @brief Compute the gradient norm of a 2D array using the Sobel filter.
 *
 * This function calculates the gradient norm (magnitude) using the Sobel filter
 * for gradient computation. Optionally, the x and y gradients can be computed
 * and stored in `p_dx` and `p_dy` respectively.
 *
 * @param  array Input 2D array.
 * @param  p_dx  Pointer to an Array where the x-gradient will be stored
 *               (optional).
 * @param  p_dy  Pointer to an Array where the y-gradient will be stored
 *               (optional).
 * @return       Array Gradient norm (magnitude).
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_norm_sobel(const Array &array,
                          Array       *p_dx = nullptr,
                          Array       *p_dy = nullptr);

/**
 * @brief Compute the gradient talus slope of a 2D array.
 *
 * The talus slope is defined as the largest elevation difference between a cell
 * and its immediate neighbors. This function computes the talus slope for the
 * given 2D array.
 *
 * @param  array Input 2D array.
 * @return       Array Gradient talus slope.
 *
 * @see          Thermal erosion: {@link thermal}.
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_talus(const Array &array);

/**
 * @brief Compute the gradient talus slope and store it in the provided array.
 *
 * This overload computes the gradient talus slope of a 2D array and stores the
 * result in the given `talus` array.
 *
 * @param array Input 2D array.
 * @param talus Output array where the gradient talus slope will be stored.
 *
 * See unit tests: @ref test_gradient.cpp
 */
void gradient_talus(const Array &array, Array &talus); ///< @overload

/**
 * @brief Compute the gradient in the x-direction of a 2D array.
 *
 * This function calculates the gradient in the x-direction (i.e., horizontal
 * gradient) for the given 2D array.
 *
 * @param  array Input 2D array.
 * @return       Array Gradient in the x-direction.
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_x(const Array &array);

/**
 * @brief Compute the gradient in the x-direction of a 2D array and store it.
 *
 * This overload calculates the gradient in the x-direction and stores the
 * result in the provided `dx` array.
 *
 * @param array Input 2D array.
 * @param dx    Output array where the gradient in the x-direction will be
 *              stored.
 */
void gradient_x(const Array &array, Array &dx); ///< @overload

/**
 * @brief Compute the gradient in the y-direction of a 2D array.
 *
 * This function calculates the gradient in the y-direction (i.e., vertical
 * gradient) for the given 2D array.
 *
 * @param  array Input 2D array.
 * @return       Array Gradient in the y-direction.
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array gradient_y(const Array &array);

/**
 * @brief Compute the gradient in the y-direction of a 2D array and store it.
 *
 * This overload calculates the gradient in the y-direction and stores the
 * result in the provided `dy` array.
 *
 * @param array Input 2D array.
 * @param dy    Output array where the gradient in the y-direction will be
 *              stored.
 */
void gradient_y(const Array &array, Array &dy); ///< @overload

/**
 * @brief Compute the Laplacian of a 2D array.
 *
 * The Laplacian is a measure of the second-order spatial derivative of the
 * input array. This function computes the Laplacian for the given 2D array.
 *
 * @param  array Input 2D array.
 * @return       Array Laplacian of the input array.
 *
 *  **Example**
 * @include ex_laplacian.cpp
 *
 * **Result**
 * @image html ex_laplacian.png
 *
 * See unit tests: @ref test_gradient.cpp
 */
Array laplacian(const Array &array);

/**
 * @brief Generates a normal map from a given 2D array.
 *
 * This function calculates the normal vectors for each element in the 2D array
 * and returns them as a tensor. The normal map is commonly used in image
 * processing and 3D graphics to represent surface orientations based on height
 * values.
 *
 * @param  array A 2D array representing the height values from which normals
 *               are computed.
 * @return       A tensor containing the normal vectors for each position in the
 *               input array.
 */
Texture normal_map(const Array &array);

/**
 * @brief Converts a normal map to a heightmap using direct summation of
 * gradients.
 *
 * This function computes a heightmap from a given normal map (`nmap`) by
 * integrating gradients derived from the normal map in two passes. The result
 * is the average of two independent integrations to reduce artifacts from
 * directional bias.
 *
 * @param  nmap A 3D tensor representing the normal map. It should have three
 *              channels (X, Y, Z) and a shape of (width, height, 3). The values
 *              are assumed to be in the range [0, 1].
 *
 * @return      An `Array` object representing the computed heightmap with the
 *              same spatial dimensions as the input normal map (width, height).
 *
 * @note
 * - The algorithm assumes the normal map values are normalized between 0 and 1.
 * - Two heightmaps (`z1` and `z2`) are computed using different traversal
 * orders, and the final result is their average to reduce directional bias.
 *
 * * **Example**
 * @include ex_normal_map_to_heightmap.cpp
 *
 * **Result**
 * @image html ex_normal_map_to_heightmap.png
 */
Array normal_map_to_heightmap(const Texture &nmap);

Array normal_map_to_heightmap_poisson(const Texture &nmap,
                                      int            iterations = 500,
                                      float          omega = 1.5f);

/**
 * @brief Solve the Poisson equation ∇²h = rhs using Gauss–Seidel iteration.
 *
 * This function assumes Dirichlet boundary conditions (height = 0 at the
 * edges). The algorithm updates interior points according to: \f[ h(x, y) = (1
 * - \omega) h(x, y)
 *           + \omega \cdot 0.25 \big(
 *               h(x+1, y) + h(x-1, y)
 *             + h(x, y+1) + h(x, y-1)
 *             - rhs(x, y)
 *             \big)
 * \f]
 *
 * @param rhs        Right-hand side (divergence of gradients).
 * @param h          Height field to update; initialized to 0 or an estimate.
 * @param iterations Number of Gauss–Seidel iterations to run.
 * @param omega      Relaxation parameter (1 = standard Gauss–Seidel, 1 < omega
 *                   < 2 = over-relaxation).
 */
void solve_poisson_gauss_seidel(const Array &rhs,
                                Array       &h,
                                int          iterations = 500,
                                float        omega = 1.0f);

/**
 * @brief Unwraps a 2D phase array to correct discontinuities in phase data.
 *
 * This function unwraps phase values in the given array, ensuring that
 * discontinuities exceeding \f$ \pi \f$ are adjusted to provide a continuous
 * phase result.
 *
 * @param  alpha A 2D array of wrapped phase values.
 * @return       A 2D array of unwrapped phase values.
 *
 * **Example**
 * @include ex_unwrap_phase.cpp
 *
 * **Result**
 * @image html ex_unwrap_phase.png
 */
Array unwrap_phase(const Array &alpha);

/* See unit tests: @ref test_gradient.cpp */
Array talus_jump_mask(const Array &z, float threshold, float sigma);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::gradient_angle_circular_smoothing */
Array gradient_angle_circular_smoothing(const Array &array,
                                        int          ir,
                                        bool         downward = false);

/**
 * @brief Compute the fractional Laplacian of an array.
 *
 * Approximates a non-local Laplacian of order @p s using a neighborhood of
 * radius @p ir. The result is unnormalized.
 *
 * @param  array Input scalar field.
 * @param  s     Fractional order (> 0 and 1 is the standard Laplacian).
 * @param  ir    Neighborhood radius.
 *
 * @return       Array Fractional Laplacian of @p array.
 *
 *  **Example**
 * @include ex_laplacian.cpp
 *
 * **Result**
 * @image html ex_laplacian.png
 */
Array laplacian_fract(const Array &array, float s, int ir);

/**
 * @brief Apply phase averaging on real and imaginary fields using a GPU kernel.
 * @param field_real Real component (input/output).
 * @param field_imag Imaginary component (input/output).
 * @param ir         Averaging radius.
 */
void phase_averaging(Array &field_real, Array &field_imag, int ir);

/**
 * @brief Compute a phase field from an input array using local orientation and
 * noise.
 *
 * @param  array            Input scalar field.
 * @param  kw               Wave vector.
 * @param  seed             Random seed.
 * @param  kp               Phase gain.
 * @param  rotate90         Rotate local angle by 90 degrees if true.
 * @param  n_kernel_samples Number of kernel samples.
 * @param  jitter           Jitter applied to sampling.
 * @param  angle_filter_ir  Radius for angle filtering.
 * @param  p_ctrl_param     Optional control parameter field.
 * @param  p_noise_x        Optional noise field (x).
 * @param  p_noise_y        Optional noise field (y).
 * @param  p_modulus        Optional output.
 * @param  bbox             Bounding box for evaluation.
 * @return                  Computed phase field.
 *
 * * **Example**
 * @include ex_phase_field.cpp
 *
 * **Result**
 * @image html ex_phase_field.png
 */
Array phase_field(const Array     &array,
                  const glm::vec2 &kw,
                  std::uint32_t    seed,
                  float            kp,
                  bool             rotate90 = false,
                  int              n_kernel_samples = 8,
                  const glm::vec2 &jitter = {0.5f, 0.5f},
                  int              angle_filter_ir = 8,
                  const Array     *p_ctrl_param = nullptr,
                  const Array     *p_noise_x = nullptr,
                  const Array     *p_noise_y = nullptr,
                  Array           *p_modulus = nullptr,
                  Array           *p_angle_jump_mask = nullptr,
                  glm::vec4        bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Compute a phase field using a global isotropic kp value.
 *
 * @param  array            Input scalar field.
 * @param  seed             Random seed.
 * @param  kp_global        Global phase gain.
 * @param  rotate90         Rotate local angle by 90 degrees if true.
 * @param  n_kernel_samples Number of kernel samples.
 * @param  jitter           Jitter applied to sampling.
 * @param  angle_filter_ir  Radius for angle filtering.
 * @param  p_ctrl_param     Optional control parameter field.
 * @param  p_noise_x        Optional noise field (x).
 * @param  p_noise_y        Optional noise field (y).
 * @param  p_modulus        Optional output.
 * @param  bbox             Bounding box for evaluation.
 * @return                  Computed phase field.
 *
 * * **Example**
 * @include ex_phase_field.cpp
 *
 * **Result**
 * @image html ex_phase_field.png
 */
Array phase_field(const Array     &array,
                  std::uint32_t    seed,
                  float            kp_global,
                  bool             rotate90 = false,
                  int              n_kernel_samples = 8,
                  const glm::vec2 &jitter = {0.5f, 0.5f},
                  int              angle_filter_ir = 8,
                  const Array     *p_ctrl_param = nullptr,
                  const Array     *p_noise_x = nullptr,
                  const Array     *p_noise_y = nullptr,
                  Array           *p_modulus = nullptr,
                  Array           *p_angle_jump_mask = nullptr,
                  glm::vec4        bbox = {0.f, 1.f, 0.f, 1.f});

Array phase_field_angle(const Array     &angle, // rads
                        const glm::vec2 &kw,
                        std::uint32_t    seed,
                        float            kp,
                        int              n_kernel_samples = 8,
                        const glm::vec2 &jitter = {0.5f, 0.5f},
                        const Array     *p_ctrl_param = nullptr,
                        const Array     *p_noise_x = nullptr,
                        const Array     *p_noise_y = nullptr,
                        Array           *p_modulus = nullptr,
                        Array           *p_angle_jump_mask = nullptr,
                        glm::vec4        bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap::gpu