/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
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

Array flow_accumulation_from_velocity_field(const Array &u,
                                            const Array &v,
                                            int          iterations)
{
  if (!validate_non_empty(u)) return Array();
  if (!validate_same_shape(u, v)) return Array();

  const glm::ivec2 shape = u.shape;
  Array            facc(shape, 1.f); // output

  const float dt = 0.5f;

  auto run = clwrapper::Run("eulerian_transport");

  run.bind_imagef("u", u.vector, shape.x, shape.y); // inputs
  run.bind_imagef("v", v.vector, shape.x, shape.y);
  run.bind_imagef("facc_in", facc.vector, shape.x, shape.y);

  run.bind_imagef("facc_out", facc.vector, shape.x, shape.y, true); // output

  run.bind_arguments(shape.x, shape.y, dt);

  for (int it = 0; it < iterations; ++it)
  {
    run.write_imagef("facc_in");
    run.execute({shape.x, shape.y});
    run.read_imagef("facc_out");
  }

  return facc;
}

} // namespace hmap::gpu
