/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/kernels.hpp"

namespace hmap::gpu
{
Array bilateral_filter(const Array &array,
                       const Array &kernel2d,
                       const Array &kernel1d,
                       float        kernel1d_value_scaling)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_non_empty(kernel2d)) return Array(array.shape);
  if (!validate_non_empty(kernel1d)) return Array(array.shape);

  const glm::ivec2 shape = array.shape;

  Array out(shape);

  auto run = clwrapper::Run("bilateral_filter");

  run.bind_imagef("array", array.vector, shape.x, shape.y);
  run.bind_imagef("kernel2d",
                  kernel2d.vector,
                  kernel2d.shape.x,
                  kernel2d.shape.y);
  run.bind_imagef("kernel1d", kernel1d.vector, kernel1d.shape.x, 1);
  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x,
                     shape.y,
                     kernel2d.shape.x,
                     kernel2d.shape.y,
                     kernel1d.shape.x,
                     kernel1d_value_scaling);

  run.write_imagef("array");
  run.execute({shape.x, shape.y});
  run.read_imagef("out");

  return out;
}

Array bilateral_filter(const Array &array, int ir, float kernel1d_value_scaling)
{
  if (!validate_non_empty(array)) return Array();

  Array kernel2d = cubic_pulse({2 * ir + 1, 2 * ir + 1});
  Array kernel1d = cubic_pulse({2 * ir + 1, 1});
  return bilateral_filter(array, kernel2d, kernel1d, kernel1d_value_scaling);
}

} // namespace hmap::gpu
