/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file array.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Declaration of the Array class for 2D floating-point arrays with
 * various mathematical operations and utilities.
 *
 * @copyright Copyright (c) 2023
 */
#pragma once

#include <functional>
#include <random>

#include <opencv2/core.hpp>

#include "highmap/algebra.hpp"
#include "highmap/colormaps.hpp"

namespace hmap
{

/**
 * @brief Array class, helper to manipulate 2D float array with "(i, j)"
 * indexing.
 *
 */
class Array
{
public:
  /**
   * @brief The shape of the array {ni, nj}.
   *
   */
  glm::ivec2 shape = glm::ivec2(0, 0);

  /**
   * @brief The underlying data storage, a vector of size shape.x * shape.y.
   *
   */
  std::vector<float> vector;

  /**
   * @brief Constructs a new Array object.
   *
   * @param shape The shape of the array {ni, nj}.
   * @param value The initial value to fill the array with.
   *
   * **Example #1**
   * @include ex_array.cpp

   * **Example #2**
   * @include ex_to_png.cpp
   *
   */
  Array();

  /// @overload
  Array(glm::ivec2 shape);

  /// @overload
  Array(glm::ivec2 shape, float value);

  /// @overload
  Array(const std::string &filename, bool flip_j = false);

  /// @overload
  Array(const std::vector<std::vector<float>> &data);

  /// @overload
  Array(const std::initializer_list<std::initializer_list<float>> &data);
  //----------------------------------------
  // overload
  //----------------------------------------

  /**
   * @brief Overloads the assignment operator for scalar assignment.
   *
   * @param  value The scalar value to assign.
   * @return       Array& Reference to the current Array object.
   */
  Array &operator=(const float value);

  /**
   * @brief Overloads the multiplication-assignment operator for scalar
   * multiplication.
   *
   * @param  value The scalar value to multiply.
   * @return       Array& Reference to the current Array object.
   */
  Array &operator*=(const float value);

  /// @overload
  Array &operator*=(const Array &array);

  /**
   * @brief Overloads the division-assignment operator for scalar division.
   *
   * @param  value The scalar value to divide by.
   * @return       Array& Reference to the current Array object.
   */
  Array &operator/=(const float value);

  /// @overload
  Array &operator/=(const Array &array);

  /**
   * @brief Overloads the addition-assignment operator for scalar addition.
   *
   * @param  value The scalar value to add.
   * @return       Array& Reference to the current Array object.
   */
  Array &operator+=(const float value);

  /// @overload
  Array &operator+=(const Array &array);

  /**
   * @brief Overloads the subtraction-assignment operator for scalar
   * subtraction.
   *
   * @param  value The scalar value to subtract.
   * @return       Array& Reference to the current Array object.
   */
  Array &operator-=(const float value);

  /// @overload
  Array &operator-=(const Array &array);

  /**
   * @brief Overloads the multiplication operator for scalar multiplication.
   *
   * @param  value The scalar value to multiply.
   * @return       Array The resulting Array after multiplication.
   */
  Array operator*(const float value) const;

  /**
   * @brief Overloads the multiplication operator for element-wise
   * multiplication with another array.
   *
   * @param  array The Array to multiply with.
   * @return       Array The resulting Array after element-wise multiplication.
   */
  Array operator*(const Array &array) const;

  /**
   * @brief Overloads the multiplication operator for scalar multiplication.
   *
   * @param  value The scalar value to multiply.
   * @param  array The Array to multiply with.
   * @return       Array The resulting Array after multiplication.
   */
  friend Array operator*(const float value, const Array &array);

  /**
   * @brief Overloads the division operator for scalar division.
   *
   * @param  value The scalar value to divide by.
   * @return       Array The resulting Array after division.
   */
  Array operator/(const float value) const;

  /**
   * @brief Overloads the division operator for element-wise division by another
   * array.
   *
   * @param  array The Array to divide by.
   * @return       Array The resulting Array after element-wise division.
   */
  Array operator/(const Array &array) const;

  /**
   * @brief Overloads the division operator for scalar division.
   *
   * @param  value The scalar value to divide by.
   * @param  array The Array to divide by.
   * @return       Array The resulting Array after division.
   */
  friend Array operator/(const float value, const Array &array);

  /**
   * @brief Overloads the addition operator for scalar addition.
   *
   * @param  value The scalar value to add.
   * @return       Array The resulting Array after addition.
   */
  Array operator+(const float value) const;

  /**
   * @brief Overloads the addition operator for element-wise addition with
   * another array.
   *
   * @param  array The Array to add.
   * @return       Array The resulting Array after element-wise addition.
   */
  Array operator+(const Array &array) const;

  /**
   * @brief Overloads the addition operator for scalar addition.
   *
   * @param  value The scalar value to add.
   * @param  array The Array to add.
   * @return       Array The resulting Array after addition.
   */
  friend Array operator+(const float value, const Array &array);

  /**
   * @brief Overloads the unary minus operator.
   *
   * @return Array The resulting Array after applying the unary minus.
   */
  Array operator-() const;

  /**
   * @brief Overloads the subtraction operator for scalar subtraction.
   *
   * @param  value The scalar value to subtract.
   * @return       Array The resulting Array after subtraction.
   */
  Array operator-(const float value) const;

  /**
   * @brief Overloads the subtraction operator for element-wise subtraction with
   * another array.
   *
   * @param  array The Array to subtract.
   * @return       Array The resulting Array after element-wise subtraction.
   */
  Array operator-(const Array &array) const;

  /**
   * @brief Overloads the subtraction operator for scalar subtraction.
   *
   * @param  value The scalar value to subtract.
   * @param  array The Array to subtract from.
   * @return       Array The resulting Array after subtraction.
   */
  friend const Array operator-(const float value, const Array &array);

  /**
   * @brief Overloads the function call operator to access the array value at
   * index (i, j).
   *
   * @param  i The row index.
   * @param  j The column index.
   * @return   float& Reference to the array value at index (i, j).
   */
  float &operator()(int i, int j)
  {
    return this->vector[j * this->shape.x + i];
  }

  /**
   * @brief Overloads the function call operator to access the array value at
   * index (i, j) (const version).
   *
   * @param  i The row index.
   * @param  j The column index.
   * @return   const float& Reference to the array value at index (i, j).
   */

  /// @overload
  const float &operator()(int i, int j) const
  {
    return this->vector[j * this->shape.x + i];
  }

  float &operator()(glm::ivec2 ij)
  {
    return this->vector[ij.y * this->shape.x + ij.x];
  }

  const float &operator()(glm::ivec2 ij) const
  {
    return this->vector[ij.y * this->shape.x + ij.x];
  }

  /**
   * @brief Overloads the function call operator to access the array value at
   * the linear index 'index' (const version).
   *
   * @param  index Linear index.
   * @return       const float& Reference to the array value.
   */

  /// @overload
  const float &operator()(int index) const
  {
    return this->vector[index];
  }

  /**
   * @brief Overloads the function call operator to access the array value at
   * the linear index 'index'.
   *
   * @param  index Linear index.
   * @return       float& Reference to the array value.
   */
  float &operator()(int index)
  {
    return this->vector[index];
  }

  //----------------------------------------
  // methods
  //----------------------------------------

  /**
   * @brief Finds the maximum value in the array and its coordinates.
   *
   * @param[out] max The maximum value found.
   * @param[out] im  The x-index of the maximum value.
   * @param[out] jm  The y-index of the maximum value.
   */
  void argmax(float &vmax, int &i, int &j) const;

  /**
   * @brief Extracts a column 'j' as a std::vector.
   *
   * @param  j The column index.
   * @return   std::vector<float> The column values as a vector.
   */
  std::vector<float> col_to_vector(int j);

  /**
   * @brief Distribute a value 'amount' around the four cells (i, j), (i + 1,
   * j), (i, j + 1), (i + 1, j + 1) by "reversing" the bilinear interpolation.
   *
   * @param i      Index.
   * @param j      Index.
   * @param u      'u' interpolation parameter, expected to be in [0, 1[.
   * @param v      'v' interpolation parameter, expected to be in [0, 1[.
   * @param amount Amount to be deposited.
   */
  void depose_amount_bilinear_at(int i, int j, float u, float v, float amount);

  void depose_amount_kernel_bilinear_at(int   i,
                                        int   j,
                                        float u,
                                        float v,
                                        int   ir,
                                        float amount);

  /**
   * @brief Distributes a specified amount of value around the cell located at
   * (i, j) using a 1D deposition kernel applied in both horizontal and vertical
   * directions.
   *
   * @param i      Index of the cell in the array where the deposition starts.
   * @param j      Index of the cell in the array where the deposition starts.
   * @param kernel 1D deposition kernel to be used for distribution. The kernel
   *               must contain an odd number of elements.
   * @param amount The amount of value to be distributed around the cell.
   */
  void depose_amount_kernel_at(int i, int j, const Array &kernel, float amount);

  /**
   * @brief Debug tool, dump some infos and generate an output file (16bits
   * grayscale)
   */
  void dump(const std::string &fname = "out.png") const;

  /**
   * @brief Debug tool, dump array content as an ASCII histogram to the consol.
   */
  void dump_histogram(int bins = 32, const std::string &msg = "") const;

  /**
   * @brief Extracts a subarray defined by the slice indices {i1, i2, j1, j2}
   *        from the original array, creating a new array. Note that i2 and j2
   * are excluded from the slice.
   *
   * @param  idx A glm::ivec4 containing the slice extent indices {i1, i2, j1,
   *             j2}.
   * @return     Array The extracted subarray.
   */
  Array extract_slice(glm::ivec4 idx) const;

  /// @overload
  Array extract_slice(int i1, int i2, int j1, int j2) const;

  /**
   * @brief Import array data from a raw binary file.
   *
   *        This function loads the array data from a binary file specified by
   * `fname`. The file should contain raw binary data that corresponds to the
   * internal representation of the array. The array's size and data type must
   * match the expected format in the file.
   *
   * @param fname The name of the file to import data from.
   */
  void from_file(const std::string &fname);

  /**
   * @brief Import array data from a numpy binary file.
   *
   * @param fname The name of the file to import data from.
   *
   * **Example**
   * @include ex_from_numpy.cpp
   */
  void from_numpy(const std::string &fname);

  /**
   * @brief Calculates the gradient in the 'x' (or 'i') direction at a specified
   * index (i, j) using a 2nd order central difference scheme.
   *
   *        This function computes the gradient in the x-direction, which
   * measures the rate of change along the x-axis. Note that this function
   * cannot be used at the borders of the data, as the central difference scheme
   * requires values from both sides of the point of interest.
   *
   * @warning The gradient cannot be computed at the borders of the data array.
   * Specifically, the function is not defined for i = 0, j = 0, i =
   * shape.x - 1, or j = shape.y - 1.
   *
   * @param  i Index along the x-direction, expected to be in the range [1,
   *           shape.x - 2].
   * @param  j Index along the y-direction, expected to be in the range [1,
   *           shape.y - 2].
   * @return   float The computed gradient value in the x-direction at the
   *           specified index.
   */
  float get_gradient_x_at(int i, int j) const;

  /**
   * @brief Calculates the gradient in the 'y' (or 'j') direction at a specified
   * index (i, j) using a 2nd order central difference scheme.
   *
   *        This function computes the gradient in the y-direction, which
   * measures the rate of change along the y-axis. Note that this function
   * cannot be used at the borders of the data, as the central difference scheme
   * requires values from both sides of the point of interest.
   *
   * @warning The gradient cannot be computed at the borders of the data array.
   * Specifically, the function is not defined for i = 0, j = 0, i =
   * shape.x - 1, or j = shape.y - 1.
   *
   * @param  i Index along the x-direction, expected to be in the range [1,
   *           shape.x - 2].
   * @param  j Index along the y-direction, expected to be in the range [1,
   *           shape.y - 2].
   * @return   float The computed gradient value in the y-direction at the
   *           specified index.
   */
  float get_gradient_y_at(int i, int j) const;

  /**
   * @brief Calculates the gradient in the 'x' (or 'i') direction at a location
   * (x, y) near the index (i, j) using bilinear interpolation.
   *
   *        This function uses bilinear interpolation to estimate the gradient
   * in the x-direction at a point within the cell defined by (i, j). The
   * gradient is based on a 2nd order central difference scheme but interpolated
   * to a finer resolution.
   *
   * @warning The gradient cannot be computed at the borders of the data array.
   * Specifically, the function is not defined for i = 0, j = 0, i =
   * shape.x - 1, or j = shape.y - 1.
   *
   * @param  i Index along the x-direction, expected to be in the range [1,
   *           shape.x - 2].
   * @param  j Index along the y-direction, expected to be in the range [1,
   *           shape.y - 2].
   * @param  u 'u' interpolation parameter, expected to be in the range [0, 1[.
   * @param  v 'v' interpolation parameter, expected to be in the range [0, 1[.
   * @return   float The interpolated gradient value in the x-direction at the
   *           specified location.
   */
  float get_gradient_x_bilinear_at(int i, int j, float u, float v) const;

  /**
   * @brief Calculates the gradient in the 'y' (or 'j') direction at a location
   * (x, y) near the index (i, j) using bilinear interpolation.
   *
   *        This function uses bilinear interpolation to estimate the gradient
   * in the y-direction at a point within the cell defined by (i, j). The
   * gradient is based on a 2nd order central difference scheme but interpolated
   * to a finer resolution.
   *
   * @warning The gradient cannot be computed at the borders of the data array.
   * Specifically, the function is not defined for i = 0, j = 0, i =
   * shape.x - 1, or j = shape.y - 1.
   *
   * @param  i Index along the x-direction, expected to be in the range [1,
   *           shape.x - 2].
   * @param  j Index along the y-direction, expected to be in the range [1,
   *           shape.y - 2].
   * @param  u 'u' interpolation parameter, expected to be in the range [0, 1[.
   * @param  v 'v' interpolation parameter, expected to be in the range [0, 1[.
   * @return   float The interpolated gradient value in the y-direction at the
   *           specified location.
   */
  float get_gradient_y_bilinear_at(int i, int j, float u, float v) const;

  /**
   * @brief Calculates the surface normal vector at the index (i, j).
   *
   *        This function computes the surface normal at the specified index (i,
   * j) in the data array. The normal vector is a three-dimensional vector (x,
   * y, z) that is perpendicular to the surface at that point.
   *
   * @param  i Index along the x-direction.
   * @param  j Index along the y-direction.
   * @return   glm::vec3 The normal vector at the specified index (i, j).
   */
  glm::vec3 get_normal_at(int i, int j) const;

  /**
   * @brief Retrieves the shape of the array.
   *
   * @return glm::ivec2 The shape {ni, nj}.
   */
  glm::ivec2 get_shape();

  /**
   * @brief Retrieves the number of bytes occupied by the array data.
   *
   *        This function returns the total size in bytes of the underlying
   * storage vector that holds the array data. This can be useful for memory
   * management or debugging purposes.
   *
   * @return size_t The size of the array data in bytes.
   */
  size_t get_sizeof() const;

  /**
   * @brief Retrieves the array value at the location (x, y) near the index (i,
   * j) using bicubic interpolation.
   *
   *        This function uses bicubic interpolation to estimate the value at a
   * point within the cell defined by (i, j). Bicubic interpolation provides a
   * smoother result compared to bilinear interpolation by taking into account
   * the values of the surrounding 4x4 grid of points, resulting in better
   * accuracy when interpolating in areas with gradients or curves.
   *
   * @warning This function cannot be used at the upper borders of the array,
   * specifically for i = shape.x - 2 or j = shape.y - 2, as bicubic
   * interpolation requires surrounding data points from a 4x4 neighborhood.
   *
   * @param  i Index along the x-direction, expected to be in the range [1,
   *           shape.x - 3].
   * @param  j Index along the y-direction, expected to be in the range [1,
   *           shape.y - 3].
   * @param  u 'u' interpolation parameter, expected to be in the range [0, 1[.
   * @param  v 'v' interpolation parameter, expected to be in the range [0, 1[.
   * @return   float The interpolated value at the specified location (x, y).
   */
  float get_value_bicubic_at(int i, int j, float u, float v) const;

  /**
   * @brief Retrieves the array value at the location (x, y) near the index (i,
   * j) using bilinear interpolation.
   *
   *        This function uses bilinear interpolation to estimate the value at a
   * point within the cell defined by (i, j). Bilinear interpolation is applied
   * based on the surrounding values to approximate the value at the specified
   * (x, y) location.
   *
   * @warning This function cannot be used at the upper borders of the array,
   * specifically for i = shape.x - 1 or j = shape.y - 1, as bilinear
   * interpolation requires surrounding data points.
   *
   * @param  i Index along the x-direction, expected to be in the range [0,
   *           shape.x - 2].
   * @param  j Index along the y-direction, expected to be in the range [0,
   *           shape.y - 2].
   * @param  u 'u' interpolation parameter, expected to be in the range [0, 1[.
   * @param  v 'v' interpolation parameter, expected to be in the range [0, 1[.
   * @return   float The interpolated value at the specified location (x, y).
   */
  float get_value_bilinear_at(int i, int j, float u, float v) const;

  /**
   * @brief Retrieves the nearest value at the location (x, y) within a bounded
   * domain.
   *
   *        This function retrieves the nearest value to the specified (x, y)
   * coordinates from the array, assuming the array covers a domain defined by
   * the bounding box `bbox`. If the coordinates (x, y) are outside the domain,
   * they are clamped to the nearest valid position within the bounding box.
   *
   * @param  x    Coordinate along the x-axis.
   * @param  y    Coordinate along the y-axis.
   * @param  bbox Bounding box defining the domain of the array. It is expected
   *              to be in the form of {xmin, xmax, ymin, ymax}.
   * @return      float The nearest value at the clamped location (x, y).
   */
  float get_value_nearest(float x, float y, glm::vec4 bbox) const;

  /**
   * @brief Retrieves the underlying data vector.
   *
   * @return std::vector<float> The vector containing the array's data.
   */
  std::vector<float> get_vector() const;

  /**
   * @brief Display various information about the array.
   *
   *        This function outputs details about the array, such as its
   * dimensions, data type, and any other relevant information. The optional
   * `msg` parameter can be used to include a custom message in the output.
   *
   * @param msg Optional message to include in the output.
   */
  void infos(const std::string &msg = "") const;

  /**
   * @brief Return the linear index corresponding to the (i, j) cell in a 2D
   * array.
   *
   *        This function calculates the linear index for a cell specified by
   * the (i, j) coordinates in a 2D array. The linear index is useful for
   * accessing elements in a flattened or 1D representation of the array.
   *
   * @param  i The 'i' index (row index).
   * @param  j The 'j' index (column index).
   * @return   int The linear index corresponding to the (i, j) cell.
   */
  int linear_index(int i, int j) const;

  /**
   * @brief Convert a linear index to its corresponding (i, j) cell coordinates.
   *
   *        This function calculates the 2D cell coordinates (i, j) in the array
   * from a given linear index `k`. This is useful for mapping between the 1D
   * and 2D representations of the array.
   *
   * @param  k The linear index.
   * @return   glm::ivec2 The (i, j) coordinates corresponding to the linear
   *           index `k`.
   */
  glm::ivec2 linear_index_reverse(int k) const;

  /**
   * @brief Return the value of the greatest element in the array.
   *
   *        This function computes the maximum value present in the array. If
   * the array is empty, the behavior of this function is undefined.
   *
   * @return float The maximum value in the array.
   */
  float max() const;

  /**
   * @brief Return the mean value of the elements in the array.
   *
   *        This function calculates the average value of all elements in the
   * array. If the array is empty, the behavior of this function is undefined.
   *
   * @return float The mean value of the elements in the array.
   */
  float mean() const;

  /**
   * @brief Return the median value of the elements in the array.
   * @return float The median value of the elements in the array.
   */
  float median() const;

  /**
   * @brief Return the value of the smallest element in the array.
   *
   *        This function computes the minimum value present in the array. If
   * the array is empty, the behavior of this function is undefined.
   *
   * @return float The minimum value in the array.
   */
  float min() const;

  /**
   * @brief Normalize array values so that the sum of all elements is equal
   * to 1.
   *
   *        This function adjusts the values in the array so that their total
   * sum is normalized to 1. This is often used in scenarios where the array
   * represents a probability distribution.
   */
  void normalize();

  /**
   * @brief Computes normalization coefficients (a, b) such that a * array + b
   * maps the values to the range [vmin, vmax].
   *
   *        This function calculates the coefficients used to normalize the
   * array values so that they fit within the specified range [vmin, vmax]. The
   * normalization is performed using the linear transformation: `a * array +
   * b`.
   *
   * @param  vmin Lower bound of the desired range.
   * @param  vmax Upper bound of the desired range.
   * @return      glm::vec2 Normalization coefficients (a, b) where `a` scales
   *              the values and `b` shifts them.
   */
  glm::vec2 normalization_coeff(float vmin = 0.f, float vmax = 1.f) const;

  /**
   * @brief Print the array values to the standard output (stdout).
   *
   *        This function prints the contents of the array in a human-readable
   * format to the standard output. This can be useful for debugging or
   * inspecting the values of the array.
   */
  void print() const;

  /**
   * @brief Return the peak-to-peak amplitude (i.e., the difference between the
   * maximum and minimum values) of the array values.
   *
   *        This function calculates the peak-to-peak amplitude of the array,
   * which is the difference between the maximum and minimum values in the
   * array.
   *
   * @return float The peak-to-peak amplitude of the array values.
   */
  float ptp() const;

  /**
   * @brief Computes the minimum and maximum values of the array.
   *
   * @return glm::vec2 A vector where x contains the minimum value and y
   *         contains the maximum value.
   */
  glm::vec2 range() const;

  /**
   * @brief Computes an approximate percentile range using a histogram.
   *
   * This function estimates the values corresponding to the given lower and
   * upper percentiles by building a histogram over the data range and
   * accumulating counts. The result is an approximation whose precision depends
   * on the number of bins.
   *
   * @param  p_low  Lower percentile in [0, 1]
   * @param  p_high Upper percentile in [0, 1]
   * @param  bins   Number of bins used for the histogram (higher = more
   *                precise)
   *
   * @return        glm::vec2 A vector where x contains the approximated lower
   *                percentile value and y contains the approximated upper
   *                percentile value
   */
  glm::vec2 range_percentile(float  p_low,
                             float  p_high,
                             size_t bins = 1024) const;

  /**
   * @brief Return an array remapped to [0, 1].
   */
  Array remapped() const;

  /**
   * @brief Return a resampled array of shape `new_shape` using bilinear
   * interpolation.
   *
   *        This function resamples the array to a new shape `new_shape` using
   * bilinear interpolation. The resampling changes the dimensions of the array
   * while preserving the data's spatial relationships.
   *
   * @param  new_shape The target shape for the resampled array.
   * @return           Array The resampled array with the specified `new_shape`.
   *
   * **Example**
   * @include ex_resample_to_shape.cpp
   *
   * **Result**
   * @image html ex_resample_to_shape.png
   */
  Array resample_to_shape(glm::ivec2 new_shape) const;

  /**
   * @brief Return a resampled array of shape `new_shape` using bicubic
   * interpolation.
   *
   *        This function resamples the array to a new shape `new_shape` using
   * bicubic interpolation. The resampling changes the dimensions of the array
   * while preserving the data's spatial relationships with a higher degree of
   * smoothness compared to bilinear interpolation.
   *
   * @param  new_shape The target shape for the resampled array.
   * @return           Array The resampled array with the specified `new_shape`.
   *
   * **Example**
   * @include ex_array_interp.cpp
   *
   * **Result**
   * @image html ex_array_interp.png
   */
  Array resample_to_shape_bicubic(glm::ivec2 new_shape) const;
  Array resample_to_shape_bilinear(glm::ivec2 new_shape) const;

  /**
   * @brief Return a resampled array of shape `new_shape` using nearest neighbor
   * interpolation.
   *
   *        This function resamples the array to a new shape `new_shape` using
   * nearest neighbor interpolation. This method selects the value from the
   * nearest neighbor in the original array, which can be faster but may produce
   * blocky artifacts compared to other interpolation methods.
   *
   * @param  new_shape The target shape for the resampled array.
   * @return           Array The resampled array with the specified `new_shape`.
   */
  Array resample_to_shape_nearest(glm::ivec2 new_shape) const;

  /**
   * @brief Return a row `i` as a `std::vector<float>`.
   *
   *        This function extracts the specified row `i` from the array and
   * returns it as a `std::vector<float>`. This is useful for operations where
   * you need to work with individual rows outside of the array structure.
   *
   * @param  i Row index to extract.
   * @return   std::vector<float> The values of the specified row as a vector.
   */
  std::vector<float> row_to_vector(int i);

  /**
   * @brief Sets a new shape for the array.
   *
   * @param new_shape The new shape of the array.
   */
  void set_shape(glm::ivec2 new_shape);

  /**
   * @brief Set the value of a slice defined by indices {i1, i2, j1, j2} to a
   * new value.
   *
   *        This function updates all elements within the specified slice of the
   * array to the given `value`. The slice is defined by the extents {i1, i2,
   * j1, j2}, where `i1` and `i2` specify the range for the rows, and `j1` and
   * `j2` specify the range for the columns.
   *
   * @param idx   Slice extent indices: {i1, i2, j1, j2}.
   * @param value The new value to set for the specified slice.
   */
  void set_slice(glm::ivec4 idx, float value);

  /// @overload
  void set_slice(glm::ivec4 idx, const Array &array);

  /**
   * @brief Return the total number of elements in the array.
   *
   *        This function returns the size of the array, which is the total
   * number of elements it contains. This can be useful for iterating over the
   * array or for checking its dimensions.
   *
   * @return int The total number of elements in the array.
   */
  int size() const;

  /**
   * @brief Return the standard deviation of all array values.
   * @return float The standard deviation of all elements in the array.
   */
  float std() const;

  /**
   * @brief Return the sum of all array values.
   *
   *        This function calculates and returns the sum of all elements in the
   * array. It can be used to obtain aggregate values or to check the overall
   * magnitude of the array's contents.
   *
   * @return float The sum of all elements in the array.
   */
  float sum() const;

  /**
   * @brief Converts a 2D `Array` to an OpenCV `cv::Mat`.
   *
   * This function converts a 2D `Array` object into an OpenCV `cv::Mat`. The
   * conversion process does not involve copying the data; instead, it uses
   * pointers to the original data in the `Array`. As a result, any
   * modifications to the `Array` will directly affect the corresponding
   * `cv::Mat` and vice versa.
   *
   * @note Since this function relies on pointers to the original data, it is
   * crucial for the user to ensure that the `Array` object remains in scope and
   * valid for the entire lifetime of the `cv::Mat`. If the `Array` is destroyed
   * or goes out of scope, the `cv::Mat` will reference invalid memory, which
   * can lead to undefined behavior.
   *
   * @return A `cv::Mat` object that shares the data with the input `Array`.
   *
   * **Example**
   * @include ex_to_cv_mat.cpp
   */
  cv::Mat to_cv_mat();

  /**
   * @brief Export the array as an OpenEXR image file.
   *
   *        This function saves the array data as an OpenEXR image file. OpenEXR
   * is a high dynamic range (HDR) image file format. The file will be created
   * or overwritten with the name specified by `fname`.
   *
   * @param fname The name of the OpenEXR file to be created or overwritten.
   *
   * **Example**
   * @include ex_to_exr.cpp
   */
  void to_exr(const std::string &fname) const;

  /**
   * @brief Export the array to a raw binary file.
   *
   *        This function writes the contents of the array to a file in raw
   * binary format. This can be useful for saving the array data to disk or for
   * transferring it between programs. The file will be created or overwritten
   * with the specified `fname`.
   *
   * @param fname The name of the file to which the array data will be written.
   */
  void to_file(const std::string &fname) const;

  /**
   * @brief Export the array to a numpy binary file.
   *
   * @param fname The name of the file to which the array data will be written.
   *
   * **Example**
   * @include ex_to_numpy.cpp
   */
  void to_numpy(const std::string &fname) const;

  /**
   * @brief Export the array as a PNG image file with a specified colormap and
   * hillshading.
   *
   *        This function saves the array data as a PNG image file. The colors
   * of the image are determined by the specified colormap (`cmap`), and
   * optional hillshading can be applied to enhance the visual representation of
   * the data. The file will be created or overwritten with the name specified
   * by `fname`.
   *
   * @param fname       The name of the PNG file to be created or overwritten.
   * @param cmap        The colormap to be used for visualizing the data. Refer
   *                    to `cmap` for available options.
   * @param hillshading Boolean flag to enable or disable hillshading. If true,
   *                    hillshading will be applied.
   *
   * **Example**
   * @include ex_to_png.cpp
   */
  void to_png(const std::string &fname,
              int                cmap,
              bool               hillshading = false,
              int                depth = CV_8U) const;

  /**
   * @brief Export the array as a grayscale PNG image file with specified bit
   * depth.
   *
   * This function saves the array data as a grayscale PNG image file. The bit
   * depth of the image is determined by the `depth` parameter, allowing for
   * either 8-bit (0-255) or 16-bit (0-65535) grayscale values. The resulting
   * image will be created or overwritten with the name specified by `fname`.
   *
   * @param fname The name of the PNG file to be created or overwritten.
   * @param depth The bit depth of the PNG image. Default is 8-bit (CV_8U). Use
   *              CV_16U for 16-bit depth if higher precision is needed.
   *
   * **Example**
   * @include ex_to_png.cpp
   */
  void to_png_grayscale(const std::string &fname,
                        int                depth = CV_8U,
                        float              vmin = 0.f,
                        float              vmax = 0.f) const;

  /**
   * @brief Export the array as a TIFF image file.
   *
   *        This function saves the array data as a TIFF (Tagged Image File
   * Format) image file. The file will be created or overwritten with the name
   * specified by `fname`.
   *
   * @param fname The name of the TIFF file to be created or overwritten.
   *
   * **Example**
   * @include ex_to_tiff.cpp
   */
  void to_tiff(const std::string &fname) const;

  /**
   * @brief Export the array as a 16-bit raw file for Unity terrain import.
   *
   *        This function saves the array data as a raw binary file with 16-bit
   * depth. The resulting file is suitable for use with Unity terrain import
   * formats. The file will be created or overwritten with the name specified by
   * `fname`.
   *
   * @param fname The name of the raw file to be created or overwritten.
   *
   * **Example**
   * @include ex_to_raw_16bit.cpp
   */
  void to_raw_16bit(const std::string &fname) const;

  /**
   * @brief Return the unique elements of the array.
   *
   *        This function extracts and returns a vector containing all unique
   * values present in the array. The values are sorted in ascending order.
   *
   * @return std::vector<float> A vector of unique values found in the array.
   */
  std::vector<float> unique_values() const;
};

// ==========================================================================
//  Functions
// ==========================================================================

/**
 * @brief Count the number of non-zero elements in the array.
 * @param  array Input array.
 * @return       Number of elements not equal to 0.
 */
size_t count_non_zero(const Array &array);

/**
 * @brief Count the number of zero elements in the array.
 * @param  array Input array.
 * @return       Number of elements equal to 0.
 */
size_t count_zero(const Array &array);

/**
 * @brief Converts an OpenCV `cv::Mat` to a 2D `Array` with optional value
 * scaling to \[0, 1\].
 *
 * This function converts an OpenCV `cv::Mat` object into a 2D `Array`. The
 * conversion process creates a new `Array` object and copies the data from the
 * `cv::Mat` into this `Array`. If the `remap_values` parameter is set to `true`
 * (the default), the values in the resulting `Array` will be scaled to the
 * interval \[0, 1\]. If set to `false`, the values will be copied directly
 * without scaling. Modifications to the `Array` will not affect the original
 * `cv::Mat` and vice versa.
 *
 * @param  mat          Reference to the OpenCV `cv::Mat` object that will be
 *                      converted.
 * @param  remap_values A boolean flag indicating whether to scale the values to
 *                      \[0, 1\]. If `true`, the values will be scaled; if
 *                      `false`, they will be copied directly. Default is
 *                      `true`.
 * @return              A 2D `Array` object containing a copy of the data from
 *                      the input
 *                      `cv::Mat`, with optional scaling to the interval \[0,
 *                      1\].
 *
 * **Example**
 * @include ex_cv_mat_to_array.cpp
 */
Array cv_mat_to_array(const cv::Mat &mat,
                      bool           remap_values = true,
                      bool           flip_j = false);

std::string get_opencv_build_information();

// ==========================================================================
// per cell operations wrapper(s)
// ==========================================================================

/// Apply a function to every cell (mutable).
template <typename Fn> inline void for_each_cell(Array &a, Fn &&fn)
{
  const int w = a.shape.x;
  const int h = a.shape.y;

  for (int j = 0; j < h; ++j)
    for (int i = 0; i < w; ++i)
      fn(i, j, a(i, j));
}

/// Apply a function to every cell (read-only).
template <typename Fn> inline void for_each_cell(const Array &a, Fn &&fn)
{
  const int w = a.shape.x;
  const int h = a.shape.y;

  for (int j = 0; j < h; ++j)
    for (int i = 0; i < w; ++i)
      fn(i, j, a(i, j));
}

} // namespace hmap
