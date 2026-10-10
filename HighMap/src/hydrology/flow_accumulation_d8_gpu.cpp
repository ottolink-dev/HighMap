/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap::gpu
{

Array flow_direction_d8(const Array &z)
{
  if (!validate_non_empty(z)) return Array();

  Array d8 = Array(z.shape);

  auto run = clwrapper::Run("flow_direction_d8");
  run.bind_buffer<float>("z", const_cast<std::vector<float> &>(z.vector));
  run.bind_buffer<float>("d8", d8.vector);
  run.bind_arguments(z.shape.x, z.shape.y);

  run.write_buffer("z");
  run.execute({z.shape.x, z.shape.y});
  run.read_buffer("d8");

  return d8;
}

} // namespace hmap::gpu
