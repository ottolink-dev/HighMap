/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap::gpu
{

// NB - A zero-sized source cannot be bound as an OpenCL image, and
// interpolating from an empty array has no defined result. Leave the target
// as-is (its constructor zero-fills it). Mirrors the CPU guard in
// interpolate_array.cpp.

// helpers

glm::vec4 helper_transform_bbox(const UvRect &uv_rect_source,
                                const UvRect &uv_rect_target)
{
  // in the OpenCL kernel, the bounding box of the source array is
  // assumed to be a unit square. Shift and rescale the target
  // bounding box according to this hypothesis
  UvRect uv_rect_target_mod = uv_rect_target;

  uv_rect_target_mod.u_min = (uv_rect_target_mod.u_min - uv_rect_source.u_min) /
                             (uv_rect_source.u_max - uv_rect_source.u_min);
  uv_rect_target_mod.u_max = (uv_rect_target_mod.u_max - uv_rect_source.u_min) /
                             (uv_rect_source.u_max - uv_rect_source.u_min);

  uv_rect_target_mod.v_min = (uv_rect_target_mod.v_min - uv_rect_source.v_min) /
                             (uv_rect_source.v_max - uv_rect_source.v_min);
  uv_rect_target_mod.v_max = (uv_rect_target_mod.v_max - uv_rect_source.v_min) /
                             (uv_rect_source.v_max - uv_rect_source.v_min);

  return uv_rect_target_mod.to_vec4();
}

// functions

void interpolate_array_bicubic(const Array &source, Array &target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  const UvRect uv_rect = UvRect::unit();

  auto run = clwrapper::Run("interpolate_array_bicubic");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(source.shape.x,
                     source.shape.y,
                     target.shape.x,
                     target.shape.y,
                     uv_rect);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_bicubic(const Array  &source,
                               Array        &target,
                               const UvRect &uv_rect_source,
                               const UvRect &uv_rect_target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  const UvRect &uv_rect_target_mod = helper_transform_bbox(uv_rect_source,
                                                           uv_rect_target);

  // compute
  auto run = clwrapper::Run("interpolate_array_bicubic");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(source.shape.x,
                     source.shape.y,
                     target.shape.x,
                     target.shape.y,
                     uv_rect_target_mod);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_bilinear(const Array &source, Array &target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  const UvRect uv_rect = UvRect::unit();

  auto run = clwrapper::Run("interpolate_array_bilinear");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, uv_rect);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_bilinear(const Array  &source,
                                Array        &target,
                                const UvRect &uv_rect_source,
                                const UvRect &uv_rect_target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  const UvRect &uv_rect_target_mod = helper_transform_bbox(uv_rect_source,
                                                           uv_rect_target);

  // compute
  auto run = clwrapper::Run("interpolate_array_bilinear");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, uv_rect_target_mod);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_lagrange(const Array &source, Array &target, int order)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  auto run = clwrapper::Run("interpolate_array_lagrange");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(source.shape.x,
                     source.shape.y,
                     target.shape.x,
                     target.shape.y,
                     order);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_nearest(const Array &source, Array &target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  const UvRect uv_rect = UvRect::unit();

  auto run = clwrapper::Run("interpolate_array_nearest");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, uv_rect);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_nearest(const Array  &source,
                               Array        &target,
                               const UvRect &uv_rect_source,
                               const UvRect &uv_rect_target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  const UvRect &uv_rect_target_mod = helper_transform_bbox(uv_rect_source,
                                                           uv_rect_target);

  // compute
  auto run = clwrapper::Run("interpolate_array_nearest");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, uv_rect_target_mod);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

} // namespace hmap::gpu
