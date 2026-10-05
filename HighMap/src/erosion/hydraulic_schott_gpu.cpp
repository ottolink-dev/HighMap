/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <array>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"

namespace hmap::gpu
{

void hydraulic_schott(Array       &z,
                      int          iterations,
                      const Array &talus,
                      float        c_erosion,
                      float        c_thermal,
                      float        c_deposition,
                      float        flow_acc_exponent,
                      float        flow_acc_exponent_depo,
                      float        flow_routing_exponent,
                      float        thermal_weight,
                      float        deposition_weight,
                      Array       *p_flow)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;
  if (p_flow && !validate_same_shape(z, *p_flow)) return;

  Array flow = p_flow ? *p_flow : Array(z.shape, 1.f);
  if (iterations <= 0)
  {
    if (p_flow) *p_flow = flow;
    return;
  }

  Array sediment(z.shape, 0.f);

  // erosion weight is always 1
  float sum_weight = 1.f + thermal_weight + deposition_weight;

  int erosion_it = (int)(10.f / sum_weight);
  int thermal_it = erosion_it + (int)(10.f * thermal_weight / sum_weight);

  auto run = clwrapper::Run("hydraulic_schott");

  glm::ivec2 shape = z.shape;
  using clwrapper::Direction;

  run.bind_imagef("z_a", z.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("flow_a", flow.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("sediment_a",
                  sediment.vector,
                  shape.x,
                  shape.y,
                  Direction::INOUT);
  run.bind_imagef("talus",
                  const_cast<std::vector<float> &>(talus.vector),
                  shape.x,
                  shape.y,
                  Direction::IN);

  run.bind_imagef("z_b", z.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("flow_b", flow.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("sediment_b",
                  sediment.vector,
                  shape.x,
                  shape.y,
                  Direction::INOUT);

  run.bind_arguments(shape.x,
                     shape.y,
                     c_erosion,
                     c_thermal,
                     c_deposition,
                     flow_acc_exponent,
                     flow_acc_exponent_depo,
                     flow_routing_exponent,
                     erosion_it,
                     thermal_it,
                     0);

  const std::array<cl::Image2D, 2> img_z = {run.get_image2d("z_a").cl_image,
                                            run.get_image2d("z_b").cl_image};
  const std::array<cl::Image2D, 2> img_flow = {
      run.get_image2d("flow_a").cl_image,
      run.get_image2d("flow_b").cl_image};
  const std::array<cl::Image2D, 2> img_sed = {
      run.get_image2d("sediment_a").cl_image,
      run.get_image2d("sediment_b").cl_image};

  int curr = 0;

  for (int it = 0; it < iterations; it++)
  {
    int next = 1 - curr;

    // inputs: z, flow, sediment
    run.set_argument(0, img_z[curr]);
    run.set_argument(1, img_flow[curr]);
    run.set_argument(2, img_sed[curr]);

    // outputs: z_new, flow_new, sediment_new
    run.set_argument(4, img_z[next]);
    run.set_argument(5, img_flow[next]);
    run.set_argument(6, img_sed[next]);

    // it index
    run.set_argument(17, it);

    run.execute_async({shape.x, shape.y});

    curr = next;
  }

  run.finish();

  run.read_imagef(curr == 0 ? "z_a" : "z_b");
  if (p_flow)
  {
    run.read_imagef(curr == 0 ? "flow_a" : "flow_b");
    *p_flow = flow;
  }
}

void hydraulic_schott(Array       &z,
                      const Array *p_mask,
                      int          iterations,
                      const Array &talus,
                      float        c_erosion,
                      float        c_thermal,
                      float        c_deposition,
                      float        flow_acc_exponent,
                      float        flow_acc_exponent_depo,
                      float        flow_routing_exponent,
                      float        thermal_weight,
                      float        deposition_weight,
                      Array       *p_flow)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_flow && !validate_same_shape(z, *p_flow)) return;

  if (!p_mask)
    hydraulic_schott(z,
                     iterations,
                     talus,
                     c_erosion,
                     c_thermal,
                     c_deposition,
                     flow_acc_exponent,
                     flow_acc_exponent_depo,
                     flow_routing_exponent,
                     thermal_weight,
                     deposition_weight,
                     p_flow);
  else
  {
    Array z_f = z;
    hydraulic_schott(z_f,
                     iterations,
                     talus,
                     c_erosion,
                     c_thermal,
                     c_deposition,
                     flow_acc_exponent,
                     flow_acc_exponent_depo,
                     flow_routing_exponent,
                     thermal_weight,
                     deposition_weight,
                     p_flow);
    z = lerp(z, z_f, *(p_mask));
  }
}

void hydraulic_schott_erosion(Array       &z,
                              int          iterations,
                              float        c_erosion,
                              float        flow_acc_exponent,
                              float        flow_routing_exponent,
                              const Array *p_moisture_map,
                              Array       *p_flow)
{
  if (!validate_non_empty(z)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;
  if (p_flow && !validate_same_shape(z, *p_flow)) return;

  Array flow = p_flow ? *p_flow : Array(z.shape, 1.f);
  if (iterations <= 0)
  {
    if (p_flow) *p_flow = flow;
    return;
  }

  Array moisture_map = p_moisture_map ? *p_moisture_map : Array(z.shape, 1.f);

  auto run = clwrapper::Run("hydraulic_schott_erosion");

  glm::ivec2 shape = z.shape;
  using clwrapper::Direction;

  run.bind_imagef("z_a", z.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("flow_a", flow.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("moisture_map",
                  moisture_map.vector,
                  shape.x,
                  shape.y,
                  Direction::IN);

  run.bind_imagef("z_b", z.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("flow_b", flow.vector, shape.x, shape.y, Direction::INOUT);

  run.bind_arguments(shape.x,
                     shape.y,
                     c_erosion,
                     flow_acc_exponent,
                     flow_routing_exponent);

  const std::array<cl::Image2D, 2> img_z = {run.get_image2d("z_a").cl_image,
                                            run.get_image2d("z_b").cl_image};
  const std::array<cl::Image2D, 2> img_flow = {
      run.get_image2d("flow_a").cl_image,
      run.get_image2d("flow_b").cl_image};

  int curr = 0;

  for (int it = 0; it < iterations; it++)
  {
    int next = 1 - curr;

    // inputs: z, flow
    run.set_argument(0, img_z[curr]);
    run.set_argument(1, img_flow[curr]);

    // outputs: z_new, flow_new
    run.set_argument(3, img_z[next]);
    run.set_argument(4, img_flow[next]);

    run.execute_async({shape.x, shape.y});

    curr = next;
  }

  run.finish();

  run.read_imagef(curr == 0 ? "z_a" : "z_b");
  if (p_flow)
  {
    run.read_imagef(curr == 0 ? "flow_a" : "flow_b");
    *p_flow = flow;
  }
}

} // namespace hmap::gpu
