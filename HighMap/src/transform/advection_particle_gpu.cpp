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
#include "highmap/interpolate/interpolate_array.hpp"
#include "highmap/math/array.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/range.hpp"
#include "highmap/texture.hpp"
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

std::vector<Array> advection_particle(
    const Array                      &z,
    const std::vector<const Array *> &advected_fields,
    int                               iterations,
    int                               nparticles,
    std::uint32_t                     seed,
    bool                              reverse,
    bool                              post_filter,
    float                             post_filter_sigma,
    float                             advection_length,
    float                             value_persistence,
    float                             inertia,
    const Array                      *p_advection_mask,
    const Array                      *p_mask)
{
  if (!validate_non_empty(z)) return {};

  for (const auto *p_field : advected_fields)
    if (!p_field || !validate_non_empty(*p_field)) return {};

  std::vector<Array> out;
  out.reserve(advected_fields.size());
  for (const auto *p_field : advected_fields)
    out.push_back(*p_field);

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
  std::vector<const Array *> ptrs;
  ptrs.reserve(advected_fields.size());
  for (const auto &f : advected_fields)
    ptrs.push_back(&f);

  return advection_particle(z,
                            ptrs,
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
                                std::vector<const Array *>{&advected_field},
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

std::vector<Array> advection_particle(
    const Array                      &z,
    const std::vector<const Array *> &advected_fields,
    int                               nparticles,
    std::uint32_t                     seed,
    bool                              reverse,
    bool                              post_filter,
    float                             post_filter_sigma,
    float                             advection_length,
    float                             value_persistence,
    float                             inertia,
    const Array                      *p_advection_mask,
    const Array                      *p_mask)
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
  std::vector<const Array *> ptrs;
  ptrs.reserve(advected_fields.size());
  for (const auto &f : advected_fields)
    ptrs.push_back(&f);

  return advection_particle(z,
                            ptrs,
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

Texture advection_particle(const Array   &z,
                           const Texture &texture,
                           int            iterations,
                           int            nparticles,
                           std::uint32_t  seed,
                           bool           reverse,
                           bool           post_filter,
                           float          post_filter_sigma,
                           float          advection_length,
                           float          value_persistence,
                           float          inertia,
                           const Array   *p_advection_mask,
                           const Array   *p_mask)
{
  if (!validate_non_empty(z)) return Texture();

  std::vector<const Array *> ptrs;
  ptrs.reserve(texture.channels.size());
  for (const auto &c : texture.channels)
    ptrs.push_back(&c);

  auto channels_out = advection_particle(z,
                                         ptrs,
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

  return Texture(channels_out);
}

Texture advection_particle(const Array   &z,
                           const Texture &texture,
                           int            nparticles,
                           std::uint32_t  seed,
                           bool           reverse,
                           bool           post_filter,
                           float          post_filter_sigma,
                           float          advection_length,
                           float          value_persistence,
                           float          inertia,
                           const Array   *p_advection_mask,
                           const Array   *p_mask)
{
  if (!validate_non_empty(z)) return Texture();

  std::vector<const Array *> ptrs;
  ptrs.reserve(texture.channels.size());
  for (const auto &c : texture.channels)
    ptrs.push_back(&c);

  auto channels_out = advection_particle(z,
                                         ptrs,
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

  return Texture(channels_out);
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
                                std::vector<const Array *>{&advected_field},
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

std::vector<Array> advection_particle(
    const Array                      &dx,
    const Array                      &dy,
    const std::vector<const Array *> &advected_fields,
    int                               nparticles,
    std::uint32_t                     seed,
    bool                              reverse,
    bool                              post_filter,
    float                             post_filter_sigma,
    float                             advection_length,
    float                             value_persistence,
    float                             inertia,
    const Array                      *p_advection_mask,
    const Array                      *p_mask)
{
  if (!validate_non_empty(dx) || !validate_same_shape(dx, dy)) return {};
  if (!validate_non_empty(advected_fields, "advected_fields")) return {};

  for (const auto *p_field : advected_fields)
  {
    if (!p_field || !validate_same_shape(dx, *p_field)) return {};
  }

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
  for (const auto *p_field : advected_fields)
  {
    concat_advected_field.insert(concat_advected_field.end(),
                                 p_field->vector.begin(),
                                 p_field->vector.end());
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
  run.write_buffer("out");
  run.write_buffer("count");

  // a zero-sized NDRange is an error on OpenCL < 2.1 (CL_INVALID_GLOBAL_WORK_
  // SIZE on macOS), e.g. with particle_density = 0: skip the launch, cells
  // with no particle hit keep their input value below
  if (nparticles > 0)
  {
    run.execute(nparticles);

    run.read_buffer("out");
    run.read_buffer("count");
  }

  // Deconcatenate results and do post-processing for each field
  std::vector<Array> out_fields;
  out_fields.reserve(num_fields);

  for (int f = 0; f < num_fields; ++f)
  {
    Array       out(shape);
    const auto &advected_field = *advected_fields[f];

    int f_offset = f * stride;

#pragma omp parallel for collapse(2)
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        int offset = f_offset + (j * shape.x + i);
        if (count(i, j))
        {
          float val = concat_out[offset] + advected_field(i, j);
          val /= (count(i, j) + 1.f);
          out(i, j) = std::clamp(val, 0.f, 1.f);
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
  std::vector<const Array *> ptrs;
  ptrs.reserve(advected_fields.size());
  for (const auto &f : advected_fields)
    ptrs.push_back(&f);

  return advection_particle(dx,
                            dy,
                            ptrs,
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

Array advection_particle_multiscale(const Array            &z,
                                    const Array            &advected_field,
                                    std::uint32_t           seed,
                                    const std::vector<int> &steps_per_level,
                                    float                   particles_ratio,
                                    float                   advection_length,
                                    float                   value_persistence,
                                    float                   inertia,
                                    float                   mix,
                                    float                   warp,
                                    bool                    reverse,
                                    bool                    post_filter,
                                    float                   post_filter_sigma,
                                    const Array            *p_advection_mask,
                                    const Array            *p_mask)
{
  if (!validate_non_empty(z)) return Array();

  auto res = advection_particle_multiscale(
      z,
      std::vector<const Array *>{&advected_field},
      seed,
      steps_per_level,
      particles_ratio,
      advection_length,
      value_persistence,
      inertia,
      mix,
      warp,
      reverse,
      post_filter,
      post_filter_sigma,
      p_advection_mask,
      p_mask);
  if (res.empty()) return Array();
  return res[0];
}

std::vector<Array> advection_particle_multiscale(
    const Array                      &z,
    const std::vector<const Array *> &advected_fields,
    std::uint32_t                     seed,
    const std::vector<int>           &steps_per_level,
    float                             particles_ratio,
    float                             advection_length,
    float                             value_persistence,
    float                             inertia,
    float                             mix,
    float                             warp,
    bool                              reverse,
    bool                              post_filter,
    float                             post_filter_sigma,
    const Array                      *p_advection_mask,
    const Array                      *p_mask)
{
  if (!validate_non_empty(z)) return {};
  if (!validate_non_empty(advected_fields, "advected_fields")) return {};

  for (const auto *p_field : advected_fields)
  {
    if (!p_field || !validate_same_shape(z, *p_field)) return {};
  }

  if (p_advection_mask && !validate_same_shape(z, *p_advection_mask)) return {};
  if (p_mask && !validate_same_shape(z, *p_mask)) return {};

  int                nlevels = static_cast<int>(steps_per_level.size());
  std::vector<Array> current_fields;
  current_fields.reserve(advected_fields.size());
  for (const auto *p_field : advected_fields)
    current_fields.push_back(*p_field);

  if (nlevels == 0) return current_fields;

  // Build halving pyramid ladder (coarsest first, final level == z.shape)
  std::vector<glm::ivec2> ladder(nlevels);
  for (int i = 0; i < nlevels; ++i)
  {
    int shift = nlevels - 1 - i;
    ladder[i] = {std::max(2, z.shape.x >> shift),
                 std::max(2, z.shape.y >> shift)};
  }

  // Iterate coarse-to-fine across resolution levels
  for (int i = 0; i < nlevels; ++i)
  {
    int level_particles = static_cast<int>(particles_ratio * ladder[i].x *
                                           ladder[i].y);
    int level_iterations = std::max(1, steps_per_level[i]);

    if (i < nlevels - 1)
    {
      // --- Coarse resolution level: compute macro-scale advection delta

      Array z_coarse = z.resample_to_shape(ladder[i]);

      std::vector<Array> fields_coarse;
      fields_coarse.reserve(current_fields.size());
      for (const auto &field : current_fields)
      {
        fields_coarse.push_back(field.resample_to_shape(ladder[i]));
      }

      Array        level_advection_mask;
      const Array *p_lvl_advection_mask = nullptr;
      if (p_advection_mask)
      {
        level_advection_mask = p_advection_mask->resample_to_shape(ladder[i]);
        p_lvl_advection_mask = &level_advection_mask;
      }

      Array        level_mask;
      const Array *p_lvl_mask = nullptr;
      if (p_mask)
      {
        level_mask = p_mask->resample_to_shape(ladder[i]);
        p_lvl_mask = &level_mask;
      }

      std::uint32_t level_seed = seed + static_cast<std::uint32_t>(i * 1000);

      std::vector<Array> fields_coarse_out = advection_particle(
          z_coarse,
          fields_coarse,
          level_iterations,
          level_particles,
          level_seed,
          reverse,
          post_filter,
          post_filter_sigma,
          advection_length,
          value_persistence,
          inertia,
          p_lvl_advection_mask,
          p_lvl_mask);

      // Compute displacement delta and upsample smoothly to full resolution
      for (size_t f = 0; f < current_fields.size(); ++f)
      {
        Array delta_coarse = fields_coarse_out[f] - fields_coarse[f];

        Array delta_full = resample_bicubic_warp(
            delta_coarse,
            z.shape,
            warp,
            seed + static_cast<std::uint32_t>(i * 37));

        current_fields[f] += std::clamp(mix, 0.f, 1.f) * delta_full;
        clamp(current_fields[f], 0.f, 1.f);
      }
    }
    else
    {
      // --- Final level: run particle advection directly on full-detail field

      std::uint32_t level_seed = seed + static_cast<std::uint32_t>(i * 1000);

      current_fields = advection_particle(z,
                                          current_fields,
                                          level_iterations,
                                          level_particles,
                                          level_seed,
                                          reverse,
                                          post_filter,
                                          post_filter_sigma,
                                          advection_length,
                                          value_persistence,
                                          inertia,
                                          p_advection_mask,
                                          p_mask);
    }
  }

  return current_fields;
}

std::vector<Array> advection_particle_multiscale(
    const Array              &z,
    const std::vector<Array> &advected_fields,
    std::uint32_t             seed,
    const std::vector<int>   &steps_per_level,
    float                     particles_ratio,
    float                     advection_length,
    float                     value_persistence,
    float                     inertia,
    float                     mix,
    float                     warp,
    bool                      reverse,
    bool                      post_filter,
    float                     post_filter_sigma,
    const Array              *p_advection_mask,
    const Array              *p_mask)
{
  std::vector<const Array *> ptrs;
  ptrs.reserve(advected_fields.size());
  for (const auto &f : advected_fields)
    ptrs.push_back(&f);

  return advection_particle_multiscale(z,
                                       ptrs,
                                       seed,
                                       steps_per_level,
                                       particles_ratio,
                                       advection_length,
                                       value_persistence,
                                       inertia,
                                       mix,
                                       warp,
                                       reverse,
                                       post_filter,
                                       post_filter_sigma,
                                       p_advection_mask,
                                       p_mask);
}

Texture advection_particle_multiscale(const Array            &z,
                                      const Texture          &texture,
                                      std::uint32_t           seed,
                                      const std::vector<int> &steps_per_level,
                                      float                   particles_ratio,
                                      float                   advection_length,
                                      float                   value_persistence,
                                      float                   inertia,
                                      float                   mix,
                                      float                   warp,
                                      bool                    reverse,
                                      bool                    post_filter,
                                      float                   post_filter_sigma,
                                      const Array            *p_advection_mask,
                                      const Array            *p_mask)
{
  if (!validate_non_empty(z)) return Texture();

  std::vector<const Array *> ptrs;
  ptrs.reserve(texture.channels.size());
  for (const auto &c : texture.channels)
    ptrs.push_back(&c);

  auto channels_out = advection_particle_multiscale(z,
                                                    ptrs,
                                                    seed,
                                                    steps_per_level,
                                                    particles_ratio,
                                                    advection_length,
                                                    value_persistence,
                                                    inertia,
                                                    mix,
                                                    warp,
                                                    reverse,
                                                    post_filter,
                                                    post_filter_sigma,
                                                    p_advection_mask,
                                                    p_mask);

  return Texture(channels_out);
}

} // namespace hmap::gpu
