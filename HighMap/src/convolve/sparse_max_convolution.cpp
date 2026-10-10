/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap::gpu
{

Array sparse_max_convolution(const Array &array,
                             const Array &kernel,
                             bool         scale_kernel_with_value,
                             float        k_smoothmax)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_non_empty(kernel)) return Array(array.shape);

  const glm::ivec2 &shape = array.shape;

  // no negative values, raises issue with atomic max in OpenCL
  const float offset = array.min();
  Array       out = array + offset;

  auto run = clwrapper::Run("sparse_max_convolution");

  run.bind_buffer<float>("array", out.vector);
  run.bind_imagef("kernel", kernel.vector, kernel.shape.x, kernel.shape.y);
  run.bind_buffer<float>("out", out.vector);

  run.write_buffer("array");
  run.write_buffer("out");

  run.bind_arguments(shape.x,
                     shape.y,
                     kernel.shape.x,
                     kernel.shape.y,
                     scale_kernel_with_value ? 1 : 0,
                     k_smoothmax);

  run.execute({shape.x, shape.y});

  run.read_buffer("out");

  return out - offset;
}

Array sparse_max_convolution(const Array &array_values,
                             const Array &array_sizes,
                             const Array &kernel,
                             float        k_smoothmax)
{
  if (!validate_non_empty(array_values)) return Array();
  if (!validate_same_shape(array_values, array_sizes))
    return Array(array_values.shape);
  if (!validate_non_empty(kernel)) return Array(array_values.shape);

  const glm::ivec2 &shape = array_values.shape;

  // no negative values, raises issue with atomic max in OpenCL
  const float offset = array_values.min();
  Array       out = array_values + offset;

  auto run = clwrapper::Run("sparse_max_convolution_with_size");

  run.bind_buffer<float>("array_values", out.vector);
  run.bind_buffer<float>("array_sizes", array_sizes.vector);
  run.bind_imagef("kernel", kernel.vector, kernel.shape.x, kernel.shape.y);
  run.bind_buffer<float>("out", out.vector);

  run.write_buffer("array_values");
  run.write_buffer("array_sizes");
  run.write_buffer("out");

  run.bind_arguments(shape.x,
                     shape.y,
                     kernel.shape.x,
                     kernel.shape.y,
                     k_smoothmax);

  run.execute({shape.x, shape.y});

  run.read_buffer("out");

  return out - offset;
}

} // namespace hmap::gpu
