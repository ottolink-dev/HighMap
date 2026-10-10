/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/transform.hpp"

namespace hmap::gpu
{

Array advection_particle(const Array  &z,
                         const Array  &advected_field,
                         int           iterations,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse,
                         bool          post_filter,
                         float         post_filter_sigma,
                         float         advection_length,
                         float         value_persistence,
                         float         inertia,
                         const Array  *p_advection_mask,
                         const Array  *p_mask)
{
  if (!validate_non_empty(z)) return Array();

  auto res = advection_particle(z,
                                std::vector<Array>{advected_field},
                                iterations,
                                nparticles,
                                seed,
                                reverse,
                                post_filter,
                                post_filter_sigma,
                                advection_length,
                                value_persistence,
                                inertia,
                                p_advection_mask,
                                p_mask);
  if (res.empty()) return Array();
  return res[0];
}

std::vector<Array> advection_particle(const Array              &z,
                                      const std::vector<Array> &advected_fields,
                                      int                       iterations,
                                      int                       nparticles,
                                      std::uint32_t             seed,
                                      bool                      reverse,
                                      bool                      post_filter,
                                      float        post_filter_sigma,
                                      float        advection_length,
                                      float        value_persistence,
                                      float        inertia,
                                      const Array *p_advection_mask,
                                      const Array *p_mask)
{
  if (!validate_non_empty(z)) return {};

  std::vector<Array> out = advected_fields;

  for (int it = 0; it < iterations; ++it)
    out = advection_particle(z,
                             out,
                             nparticles,
                             seed,
                             reverse,
                             post_filter,
                             post_filter_sigma,
                             advection_length,
                             value_persistence,
                             inertia,
                             p_advection_mask,
                             p_mask);

  return out;
}

Array advection_particle(const Array  &z,
                         const Array  &advected_field,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse,
                         bool          post_filter,
                         float         post_filter_sigma,
                         float         advection_length,
                         float         value_persistence,
                         float         inertia,
                         const Array  *p_advection_mask,
                         const Array  *p_mask)
{
  if (!validate_non_empty(z)) return Array();

  auto res = advection_particle(z,
                                std::vector<Array>{advected_field},
                                nparticles,
                                seed,
                                reverse,
                                post_filter,
                                post_filter_sigma,
                                advection_length,
                                value_persistence,
                                inertia,
                                p_advection_mask,
                                p_mask);
  if (res.empty()) return Array();
  return res[0];
}

std::vector<Array> advection_particle(const Array              &z,
                                      const std::vector<Array> &advected_fields,
                                      int                       nparticles,
                                      std::uint32_t             seed,
                                      bool                      reverse,
                                      bool                      post_filter,
                                      float        post_filter_sigma,
                                      float        advection_length,
                                      float        value_persistence,
                                      float        inertia,
                                      const Array *p_advection_mask,
                                      const Array *p_mask)
{
  if (!validate_non_empty(z)) return {};

  Array dx = -hmap::gradient_x(z);
  Array dy = -hmap::gradient_y(z);

  return advection_particle(dx,
                            dy,
                            advected_fields,
                            nparticles,
                            seed,
                            reverse,
                            post_filter,
                            post_filter_sigma,
                            advection_length,
                            value_persistence,
                            inertia,
                            p_advection_mask,
                            p_mask);
}

Array advection_particle(const Array  &dx,
                         const Array  &dy,
                         const Array  &advected_field,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse,
                         bool          post_filter,
                         float         post_filter_sigma,
                         float         advection_length,
                         float         value_persistence,
                         float         inertia,
                         const Array  *p_advection_mask,
                         const Array  *p_mask)
{
  if (!validate_non_empty(dx)) return Array();

  auto res = advection_particle(dx,
                                dy,
                                std::vector<Array>{advected_field},
                                nparticles,
                                seed,
                                reverse,
                                post_filter,
                                post_filter_sigma,
                                advection_length,
                                value_persistence,
                                inertia,
                                p_advection_mask,
                                p_mask);
  if (res.empty()) return Array();
  return res[0];
}

std::vector<Array> advection_particle(const Array              &dx,
                                      const Array              &dy,
                                      const std::vector<Array> &advected_fields,
                                      int                       nparticles,
                                      std::uint32_t             seed,
                                      bool                      reverse,
                                      bool                      post_filter,
                                      float        post_filter_sigma,
                                      float        advection_length,
                                      float        value_persistence,
                                      float        inertia,
                                      const Array *p_advection_mask,
                                      const Array *p_mask)
{
  if (!validate_non_empty(dx) || !validate_same_shape(dx, dy)) return {};
  if (!validate_non_empty(advected_fields, "advected_fields")) return {};

  for (const auto &field : advected_fields)
    if (!validate_same_shape(dx, field)) return {};

  if (p_advection_mask && !validate_same_shape(dx, *p_advection_mask))
    return {};
  if (p_mask && !validate_same_shape(dx, *p_mask)) return {};

  auto run = clwrapper::Run("advection_particle");

  glm::ivec2 shape = dx.shape;
  int        num_fields = static_cast<int>(advected_fields.size());
  int        stride = shape.x * shape.y;

  // Concatenate input fields
  std::vector<float> concat_advected_field;
  concat_advected_field.reserve(num_fields * stride);
  for (const auto &field : advected_fields)
  {
    concat_advected_field.insert(concat_advected_field.end(),
                                 field.vector.begin(),
                                 field.vector.end());
  }

  std::vector<float> concat_out(num_fields * stride, 0.f);
  Array              count(shape);

  run.bind_buffer<float>("advected_field", concat_advected_field);
  run.bind_buffer<float>("dx", dx.vector);
  run.bind_buffer<float>("dy", dy.vector);
  run.bind_buffer<float>("out", concat_out);
  run.bind_buffer<float>("count", count.vector);
  helper_bind_optional_buffer(run, "advection_mask", p_advection_mask);

  run.bind_arguments(shape.x,
                     shape.y,
                     nparticles,
                     seed,
                     reverse ? -1.f : 1.f,
                     advection_length,
                     value_persistence,
                     inertia,
                     p_advection_mask ? 1 : 0,
                     num_fields);

  run.write_buffer("advected_field");
  run.write_buffer("dx");
  run.write_buffer("dy");

  run.execute(nparticles);

  run.read_buffer("out");
  run.read_buffer("count");

  // Deconcatenate results and do post-processing for each field
  std::vector<Array> out_fields;
  out_fields.reserve(num_fields);

  for (int f = 0; f < num_fields; ++f)
  {
    Array       out(shape);
    const auto &advected_field = advected_fields[f];

    for (int j = 0; j < shape.x; ++j)
      for (int i = 0; i < shape.y; ++i)
      {
        int offset = f * stride + (j * shape.x + i);
        if (count(i, j))
        {
          out(i, j) = concat_out[offset] + advected_field(i, j);
          out(i, j) /= (count(i, j) + 1.f);
          out(i, j) = std::clamp(out(i, j), 0.f, 1.f);
        }
        else
        {
          out(i, j) = advected_field(i, j);
        }
      }

    // post-processing
    if (post_filter)
    {
      int max_it = 1;
      hmap::laplace(out, p_advection_mask, post_filter_sigma, max_it);
    }

    if (p_mask)
      out_fields.push_back(lerp(out, advected_field, *p_mask));
    else
      out_fields.push_back(out);
  }

  return out_fields;
}

} // namespace hmap::gpu
