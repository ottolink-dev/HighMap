/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <utility>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/range.hpp"

#include <source_location>

namespace hmap
{

void add_kernel(Array &array, const Array &kernel, int ic, int jc)
{
  if (!validate_non_empty(array) || !validate_non_empty(kernel)) return;

  // truncate kernel to make it fit into the heightmap array
  int nk_i0 = (int)(std::floor(0.5f * kernel.shape.x)); // left
  int nk_i1 = kernel.shape.x - nk_i0;                   // right
  int nk_j0 = (int)(std::floor(0.5f * kernel.shape.y));
  int nk_j1 = kernel.shape.y - nk_j0;

  int ik0 = std::max(0, nk_i0 - ic);
  int jk0 = std::max(0, nk_j0 - jc);
  int ik1 = std::min(kernel.shape.x,
                     kernel.shape.x - (ic + nk_i1 - array.shape.x));
  int jk1 = std::min(kernel.shape.y,
                     kernel.shape.y - (jc + nk_j1 - array.shape.y));

  // where it goes in the array
  int i0 = std::max(ic - nk_i0, 0);
  int j0 = std::max(jc - nk_j0, 0);

  for (int j = jk0; j < jk1; j++)
    for (int i = ik0; i < ik1; i++)
      array(i - ik0 + i0, j - jk0 + j0) += kernel(i, j);
}

void add_kernel_maximum_smooth(Array       &array,
                               const Array &kernel,
                               float        k_smooth,
                               int          ic,
                               int          jc)
{
  if (!validate_non_empty(array) || !validate_non_empty(kernel)) return;

  // truncate kernel to make it fit into the heightmap array
  int nk_i0 = (int)(std::floor(0.5f * kernel.shape.x)); // left
  int nk_i1 = kernel.shape.x - nk_i0;                   // right
  int nk_j0 = (int)(std::floor(0.5f * kernel.shape.y));
  int nk_j1 = kernel.shape.y - nk_j0;

  int ik0 = std::max(0, nk_i0 - ic);
  int jk0 = std::max(0, nk_j0 - jc);
  int ik1 = std::min(kernel.shape.x,
                     kernel.shape.x - (ic + nk_i1 - array.shape.x));
  int jk1 = std::min(kernel.shape.y,
                     kernel.shape.y - (jc + nk_j1 - array.shape.y));

  // where it goes in the array
  int i0 = std::max(ic - nk_i0, 0);
  int j0 = std::max(jc - nk_j0, 0);

  for (int j = jk0; j < jk1; j++)
    for (int i = ik0; i < ik1; i++)
    {
      float v_prev = array(i - ik0 + i0, j - jk0 + j0);
      float v_new = kernel(i, j);
      array(i - ik0 + i0,
            j - jk0 + j0) = maximum_smooth(v_prev, v_new, k_smooth);
    }
}

Array hstack(const Array &array1, const Array &array2) // friend function
{
  if (!validate_non_empty(array1) || !validate_non_empty(array2))
    return Array();

  if (array1.shape.y != array2.shape.y)
  {
    hmap::log::warn(std::source_location::current(),
                    "hstack height mismatch: lhs is ({}, {}), rhs is ({}, {})",
                    array1.shape.x,
                    array1.shape.y,
                    array2.shape.x,
                    array2.shape.y);
    return Array();
  }

  Array array_out = Array(
      glm::ivec2(array1.shape.x + array2.shape.x, array1.shape.y));

  for (int j = 0; j < array1.shape.y; j++)
    for (int i = 0; i < array1.shape.x; i++)
      array_out(i, j) = array1(i, j);

  for (int j = 0; j < array1.shape.y; j++)
    for (int i = 0; i < array2.shape.x; i++)
      array_out(i + array1.shape.x, j) = array2(i, j);

  return array_out;
}

void swap(Array &a, Array &b)
{
  if (!validate_non_empty(a) || !validate_same_shape(a, b)) return;

  for (int j = 0; j < a.shape.y; j++)
    for (int i = 0; i < a.shape.x; i++)
      std::swap(a(i, j), b(i, j));
}

Array vstack(const Array &array1, const Array &array2) // friend function
{
  if (!validate_non_empty(array1) || !validate_non_empty(array2))
    return Array();

  if (array1.shape.x != array2.shape.x)
  {
    hmap::log::warn(std::source_location::current(),
                    "vstack width mismatch: lhs is ({}, {}), rhs is ({}, {})",
                    array1.shape.x,
                    array1.shape.y,
                    array2.shape.x,
                    array2.shape.y);
    return Array();
  }

  Array array_out = Array(
      glm::ivec2(array1.shape.x, array1.shape.y + array2.shape.y));

  for (int i = 0; i < array1.shape.x; i++)
  {
    for (int j = 0; j < array1.shape.y; j++)
      array_out(i, j) = array1(i, j);

    for (int j = 0; j < array2.shape.y; j++)
      array_out(i, j + array1.shape.y) = array2(i, j);
  }

  return array_out;
}

} // namespace hmap
