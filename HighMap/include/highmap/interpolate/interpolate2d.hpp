/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file interpolate2d.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for 2D interpolation methods.
 *
 * This file provides declarations for functions and enumerations related to 2D
 * interpolation. It supports different interpolation methods, including
 * Delaunay triangulation and nearest point interpolation. These functions allow
 * for the interpolation of values on a 2D grid using various techniques.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */

#pragma once
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"

extern "C"
{
  struct delaunay;
  struct nnai;
}

namespace hmap
{

/**
 * @enum InterpolationMethod2D
 * @brief Enumeration of 2D interpolation methods.
 *
 * This enum defines the available methods for 2D interpolation, such as
 * Delaunay triangulation and nearest point interpolation.
 */
enum InterpolationMethod2D : int
{
  ITP2D_DELAUNAY, ///< Delaunay triangulation method for 2D interpolation.
  ITP2D_NEAREST,  ///< Nearest point method for 2D interpolation.
  ITP2D_IDW,      ///< Inverse Distance Weighting.
  ITP2D_GAUSSIAN, ///< Gaussian Distance Weighting.
  ITP2D_NNI,      ///< Natural Neighbor Interpolation.
  ITP2D_DELAUNAY_GRADIENT, ///< Delaunay triangulation + linear gradient.
};

class NaturalNeighborInterpolator
{
public:
  NaturalNeighborInterpolator() = default;
  ~NaturalNeighborInterpolator();

  NaturalNeighborInterpolator(const NaturalNeighborInterpolator &) = delete;
  NaturalNeighborInterpolator &operator=(const NaturalNeighborInterpolator &) =
      delete;

  NaturalNeighborInterpolator(NaturalNeighborInterpolator &&other) noexcept;
  NaturalNeighborInterpolator &operator=(
      NaturalNeighborInterpolator &&other) noexcept;

  void build(const std::vector<float> &xin, const std::vector<float> &yin);

  void setup_output_points(const std::vector<float> &x,
                           const std::vector<float> &y);

  void setup_output_points(const std::vector<double> &x,
                           const std::vector<double> &y);

  void setup_output_points(std::vector<double> &&x, std::vector<double> &&y);

  void interpolate(const std::vector<float> &values_in,
                   std::vector<float>       &values_out) const;

private:
  nnai               *handle = nullptr;
  delaunay           *d = nullptr;
  std::vector<double> xout;
  std::vector<double> yout;
  size_t              nout = 0;
  size_t              nin = 0;
};

/**
 * @brief Compute the bilinear interpolated value from four input values.
 *
 * This function calculates the interpolated value at a point within a grid
 * using bilinear interpolation based on four surrounding values.
 *
 * @param  f00 Value at (u, v) = (0, 0).
 * @param  f10 Value at (u, v) = (1, 0).
 * @param  f01 Value at (u, v) = (0, 1).
 * @param  f11 Value at (u, v) = (1, 1).
 * @param  u   The interpolation parameter in the x-direction, expected in [0,
 *             1).
 * @param  v   The interpolation parameter in the y-direction, expected in [0,
 *             1).
 * @return     float The bilinear interpolated value.
 */
[[nodiscard]] inline constexpr float bilinear_interp(float f00,
                                                     float f10,
                                                     float f01,
                                                     float f11,
                                                     float u,
                                                     float v) noexcept
{
  float a10 = f10 - f00;
  float a01 = f01 - f00;
  float a11 = f11 - f10 - f01 + f00;
  return f00 + a10 * u + a01 * v + a11 * u * v;
}

[[nodiscard]] inline constexpr float cubic_interpolate(const float p[4],
                                                       float       x) noexcept
{
  return p[1] + 0.5f * x *
                    (p[2] - p[0] +
                     x * (2.f * p[0] - 5.f * p[1] + 4.f * p[2] - p[3] +
                          x * (3.f * (p[1] - p[2]) + p[3] - p[0])));
}

/**
 * @brief Perform harmonic interpolation on a 2D array using the Successive
 * Over-Relaxation (SOR) method.
 *
 * This function solves the discrete Laplace equation on a regular grid by
 * iteratively updating the values of an input array. Values marked as fixed in
 * the @p mask_fixed_values remain unchanged throughout the process. The
 * algorithm stops when either the maximum number of iterations is reached or
 * the change between successive iterations falls below the specified tolerance.
 *
 * @param  array             Input 2D array providing the initial guess for the
 *                           solution.
 * @param  mask_fixed_values Mask of the same shape as @p array. Cells with a
 *                           value greater than zero indicate fixed points that
 *                           must remain unchanged during interpolation.
 * @param  iterations_max    Maximum number of SOR iterations to perform.
 * @param  tolerance         Convergence criterion: the algorithm stops if the
 *                           maximum absolute update between iterations is less
 *                           than this value.
 * @param  omega             Relaxation factor (1 < omega < 2 recommended).
 *                           Values closer to 2 accelerate convergence, but
 *                           overly large values may cause divergence.
 *
 * @return                   A new 2D array containing the interpolated
 *                           solution.
 *
 * @note
 *  - The algorithm updates only the interior points (indices `[1..nx-2,
 * 1..ny-2]`).
 *  - Cells where `mask_fixed_values(i, j) > 0` remain unchanged.
 *  - A relaxation factor of 1 corresponds to the standard Jacobi/Gauss-Seidel
 * update.
 */
Array harmonic_interpolation(const Array &array,
                             const Array &mask_fixed_values,
                             int          iterations_max = 500,
                             float        tolerance = 1e-5f,
                             float        omega = 0.f);

/**
 * @brief Perform harmonic interpolation on a 2D array with 2D variable
 * diffusion coefficients (Dx, Dy) using the Successive Over-Relaxation (SOR)
 * method.
 *
 * Solves the steady-state anisotropic diffusion equation:
 * \f[
 * \nabla \cdot (\mathbf{D} \nabla u) = \frac{\partial}{\partial x}\left(D_x
 * \frac{\partial u}{\partial x}\right) +
 * \frac{\partial}{\partial y}\left(D_y \frac{\partial u}{\partial y}\right) = 0
 * \f]
 *
 * @param  array             Input 2D array providing the initial guess.
 * @param  mask_fixed_values Mask where values > 0 indicate fixed boundary
 *                           constraints.
 * @param  dx                Array of diffusion coefficients in the x-direction.
 * @param  dy                Array of diffusion coefficients in the y-direction.
 * @param  iterations_max    Maximum number of SOR iterations.
 * @param  tolerance         Convergence tolerance.
 * @param  omega             Relaxation factor (0 = auto-computed optimal).
 * @return                   Interpolated array.
 */
Array harmonic_interpolation(const Array &array,
                             const Array &mask_fixed_values,
                             const Array &dx,
                             const Array &dy,
                             int          iterations_max = 500,
                             float        tolerance = 1e-5f,
                             float        omega = 0.f);

/**
 * @brief Generic 2D interpolation function.
 *
 * This function performs interpolation on a 2D grid using the specified
 * interpolation method. It can optionally apply noise and stretching to the
 * input data before interpolation.
 *
 * @param  shape                Output array shape.
 * @param  x                    x coordinates of the input values.
 * @param  y                    y coordinates of the input values.
 * @param  values               Input values at (x, y).
 * @param  interpolation_method Interpolation method (see @ref
 *                              InterpolationMethod2D).
 * @param  p_noise_x            Pointer to the input noise array in the x
 *                              direction (optional).
 * @param  p_noise_y            Pointer to the input noise array in the y
 *                              direction (optional).
 * @param  bbox                 Domain bounding box (default: {0.f, 1.f, 0.f,
 *                              1.f}).
 * @return                      Array Output array with interpolated values.
 *
 * **Example**
 * @include ex_interpolate2d.cpp
 *
 * **Result**
 * @image html ex_interpolate2d.png
 */
Array interpolate2d(glm::ivec2                shape,
                    const std::vector<float> &x,
                    const std::vector<float> &y,
                    const std::vector<float> &values,
                    InterpolationMethod2D     interpolation_method,
                    const Array              *p_noise_x = nullptr,
                    const Array              *p_noise_y = nullptr,
                    glm::vec4                 bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief 2D interpolation using the Delaunay triangulation method.
 *
 * This function performs 2D interpolation by generating a Delaunay
 * triangulation of the input points and interpolating the values within each
 * triangle.
 *
 * @param  shape     Output array shape.
 * @param  x         x coordinates of the input values.
 * @param  y         y coordinates of the input values.
 * @param  values    Input values at (x, y).
 * @param  p_noise_x Pointer to the input noise array in the x direction
 *                   (optional).
 * @param  p_noise_y Pointer to the input noise array in the y direction
 *                   (optional).
 * @param  bbox      Domain bounding box (default: {0.f, 1.f, 0.f, 1.f}).
 * @return           Array Output array with interpolated values.
 */
Array interpolate2d_delaunay(glm::ivec2                shape,
                             const std::vector<float> &x,
                             const std::vector<float> &y,
                             const std::vector<float> &values,
                             const Array              *p_noise_x = nullptr,
                             const Array              *p_noise_y = nullptr,
                             glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
                             float     fill_value = 0.f);

/**
 * @brief 2D interpolation using a smoother version of the Delaunay
 * triangulation method.
 */
Array interpolate2d_delaunay_gradient(glm::ivec2                shape,
                                      const std::vector<float> &x,
                                      const std::vector<float> &y,
                                      const std::vector<float> &values,
                                      const Array *p_noise_x = nullptr,
                                      const Array *p_noise_y = nullptr,
                                      glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                                      float        fill_value = 0.f,
                                      float        gradient_scaling = 1.f);

/**
 * @brief 2D interpolation using the Gaussian kernel method.
 */
Array interpolate2d_gaussian(glm::ivec2                shape,
                             const std::vector<float> &x,
                             const std::vector<float> &y,
                             const std::vector<float> &values,
                             const Array              *p_noise_x = nullptr,
                             const Array              *p_noise_y = nullptr,
                             glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
                             float     sigma = 0.05f,
                             float     radius = 0.f);

/**
 * @brief 2D interpolation using the IDW method.
 */
Array interpolate2d_idw(glm::ivec2                shape,
                        const std::vector<float> &x,
                        const std::vector<float> &y,
                        const std::vector<float> &values,
                        const Array              *p_noise_x = nullptr,
                        const Array              *p_noise_y = nullptr,
                        glm::vec4                 bbox = {0.f, 1.f, 0.f, 1.f},
                        float                     distance_exp = 2.f,
                        float                     radius = 0.f);

/**
 * @brief 2D interpolation using the nearest neighbor method.
 *
 * This function performs 2D interpolation by assigning the value of the nearest
 * point to each point in the output grid.
 *
 * @param  shape     Output array shape.
 * @param  x         x coordinates of the input values.
 * @param  y         y coordinates of the input values.
 * @param  values    Input values at (x, y).
 * @param  p_noise_x Pointer to the input noise array in the x direction
 *                   (optional).
 * @param  p_noise_y Pointer to the input noise array in the y direction
 *                   (optional).
 * @param  bbox      Domain bounding box (default: {0.f, 1.f, 0.f, 1.f}).
 * @return           Array Output array with interpolated values.
 */
Array interpolate2d_nearest(glm::ivec2                shape,
                            const std::vector<float> &x,
                            const std::vector<float> &y,
                            const std::vector<float> &values,
                            const Array              *p_noise_x = nullptr,
                            const Array              *p_noise_y = nullptr,
                            glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief 2D interpolation using the Natural Neighbor Interpolation method.
 */
Array interpolate2d_nni(glm::ivec2                shape,
                        const std::vector<float> &x,
                        const std::vector<float> &y,
                        const std::vector<float> &values,
                        const Array              *p_noise_x = nullptr,
                        const Array              *p_noise_y = nullptr,
                        glm::vec4                 bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap

namespace hmap::gpu
{

/*! @brief Perform harmonic interpolation on GPU using Red-Black Successive
   Over-Relaxation (RB-SOR).
 *
 * @param  array             Input 2D array providing the initial guess.
 * @param  mask_fixed_values Mask array where values > 0 are fixed boundary
 *                           constraints.
 * @param  iterations_max    Maximum number of iterations.
 * @param  tolerance         Convergence tolerance for early exit (default:
 *                           1e-5). If the maximum update in an iteration check
 * falls below this threshold, the algorithm stops.
 * @param  omega             Relaxation parameter (1.0 < omega < 2.0). If set to
 *                           0.0f (default), the optimal analytical relaxation
 *                           factor for the grid shape is automatically
 *                           computed.
 * @return                   Interpolated array.
 */
Array harmonic_interpolation(const Array &array,
                             const Array &mask_fixed_values,
                             int          iterations_max = 500,
                             float        tolerance = 1e-5f,
                             float        omega = 0.f);

/*! @brief Perform harmonic interpolation on GPU with 2D variable diffusion
   coefficients (Dx, Dy) using Red-Black SOR.
 *
 * @param  array             Input 2D array providing the initial guess.
 * @param  mask_fixed_values Mask where values > 0 are fixed boundary
 *                           constraints.
 * @param  dx                Array of diffusion coefficients in x-direction.
 * @param  dy                Array of diffusion coefficients in y-direction.
 * @param  iterations_max    Maximum number of iterations.
 * @param  tolerance         Convergence tolerance.
 * @param  omega             Relaxation factor (0 = auto-computed optimal).
 * @return                   Interpolated array.
 */
Array harmonic_interpolation(const Array &array,
                             const Array &mask_fixed_values,
                             const Array &dx,
                             const Array &dy,
                             int          iterations_max = 500,
                             float        tolerance = 1e-5f,
                             float        omega = 0.f);

} // namespace hmap::gpu