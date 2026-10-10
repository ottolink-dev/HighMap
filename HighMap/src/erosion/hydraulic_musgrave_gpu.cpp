/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/erosion/deprecated.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap::gpu
{

//----------------------------------------------------------------------
// Main operator
//----------------------------------------------------------------------

void hydraulic_musgrave(Array &z,
                        Array &moisture_map,
                        int    iterations,
                        float  c_capacity,
                        float  c_erosion,
                        float  c_deposition,
                        float  water_level,
                        float  evap_rate)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, moisture_map)) return;

  Array s = constant(z.shape, 0.f);     // sediment level
  Array w = water_level * moisture_map; // initial moisture map

  Array z_new = z;
  Array w_new = w;
  Array s_new = s;

  auto run = clwrapper::Run("hydraulic_musgrave");

  glm::ivec2 shape = z.shape;

  run.bind_imagef("z", z.vector, shape.x, shape.y);
  run.bind_imagef("w", w.vector, shape.x, shape.y);
  run.bind_imagef("s", s.vector, shape.x, shape.y);
  run.bind_imagef("moisture_map", moisture_map.vector, shape.x, shape.y);

  run.bind_imagef("z_new", z_new.vector, shape.x, shape.y, true);
  run.bind_imagef("w_new", w_new.vector, shape.x, shape.y, true);
  run.bind_imagef("s_new", s_new.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x,
                     shape.y,
                     c_capacity,
                     c_erosion,
                     c_deposition,
                     water_level,
                     evap_rate);

  for (int it = 0; it < iterations; it++)
  {
    run.execute({shape.x, shape.y});

    // read outputs from device
    run.read_imagef("z_new");
    run.read_imagef("w_new");
    run.read_imagef("s_new");

    // copy to inputs for next iteration
    z = z_new;
    w = w_new;
    s = s_new;

    run.write_imagef("z");
    run.write_imagef("w");
    run.write_imagef("s");
  }
}

void hydraulic_musgrave(Array &z,
                        int    iterations,
                        float  c_capacity,
                        float  c_erosion,
                        float  c_deposition,
                        float  water_level,
                        float  evap_rate)
{
  if (!validate_non_empty(z)) return;

  Array moisture_map = constant(z.shape, 1.f);
  hmap::gpu::hydraulic_musgrave(z,
                                moisture_map,
                                iterations,
                                c_capacity,
                                c_erosion,
                                c_deposition,
                                water_level,
                                evap_rate);
}

} // namespace hmap::gpu
