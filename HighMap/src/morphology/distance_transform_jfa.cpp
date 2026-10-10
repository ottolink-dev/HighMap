/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/operator.hpp"

namespace hmap::gpu
{

Array distance_transform_jfa(const Array &array, bool return_squared_distance)
{
  if (!validate_non_empty(array)) return Array();

  const glm::ivec2 shape = array.shape;

  // --- prepare

  Array edt(shape);

  // define seed indices
  Array i_prev(shape, -1);
  Array j_prev(shape, -1);
  Array i_next(shape);
  Array j_next(shape);

  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
      if (array(i, j) != 0.f)
      {
        // foreground is own seed index (for background it is (-1,
        // -1) for now)
        i_prev(i, j) = i;
        j_prev(i, j) = j;
      }

  // --- jump flooding

  // initial step = max(width, height) / 2, iterate halving until step >= 1
  int max_dim = std::max(shape.x, shape.y);
  int step = 1;

  // find largest power-of-two-like jump: integer division by 2 until 0
  {
    step = 1;
    while (step * 2 < max_dim)
      step *= 2;
  }

  // jfa kernel
  auto run = clwrapper::Run("jump_flooding");

  run.bind_buffer<float>("i_prev", i_prev.vector);
  run.bind_buffer<float>("j_prev", j_prev.vector);
  run.bind_buffer<float>("i_next", i_next.vector);
  run.bind_buffer<float>("j_next", j_next.vector);
  run.bind_buffer<float>("edt", edt.vector);
  run.bind_arguments(shape.x, shape.y, step);

  while (step > 0)
  {
    run.write_buffer("i_prev");
    run.write_buffer("j_prev");

    run.set_argument(7, step);

    run.execute({shape.x, shape.y});

    run.read_buffer("i_next");
    run.read_buffer("j_next");

    // swap buffers
    swap(i_prev, i_next);
    swap(j_prev, j_next);

    // next it
    step /= 2;
  }

  // --- output

  run.read_buffer("edt");

  if (return_squared_distance)
    return edt;
  else
    return sqrt(edt);
}

} // namespace hmap::gpu
