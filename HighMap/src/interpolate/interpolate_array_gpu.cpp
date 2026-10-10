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

glm::vec4 helper_transform_bbox(const glm::vec4 &bbox_source,
                                const glm::vec4 &bbox_target)
{
  // in the OpenCL kernel, the bounding box of the source array is
  // assumed to be a unit square. Shift and rescale the target
  // bounding box according to this hypothesis
  glm::vec4 bbox_target_mod = bbox_target;

  bbox_target_mod.x = (bbox_target_mod.x - bbox_source.x) /
                      (bbox_source.y - bbox_source.x);
  bbox_target_mod.y = (bbox_target_mod.y - bbox_source.x) /
                      (bbox_source.y - bbox_source.x);

  bbox_target_mod.z = (bbox_target_mod.z - bbox_source.z) /
                      (bbox_source.w - bbox_source.z);
  bbox_target_mod.w = (bbox_target_mod.w - bbox_source.z) /
                      (bbox_source.w - bbox_source.z);

  return bbox_target_mod;
}

// functions

void interpolate_array_bicubic(const Array &source, Array &target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  glm::vec4 bbox(0.f, 1.f, 0.f, 1.f);

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
                     bbox);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_bicubic(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  glm::vec4 bbox_target_mod = helper_transform_bbox(bbox_source, bbox_target);

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
                     bbox_target_mod);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_bilinear(const Array &source, Array &target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  glm::vec4 bbox(0.f, 1.f, 0.f, 1.f);

  auto run = clwrapper::Run("interpolate_array_bilinear");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, bbox);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_bilinear(const Array     &source,
                                Array           &target,
                                const glm::vec4 &bbox_source,
                                const glm::vec4 &bbox_target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  glm::vec4 bbox_target_mod = helper_transform_bbox(bbox_source, bbox_target);

  // compute
  auto run = clwrapper::Run("interpolate_array_bilinear");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, bbox_target_mod);

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

  glm::vec4 bbox(0.f, 1.f, 0.f, 1.f);

  auto run = clwrapper::Run("interpolate_array_nearest");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, bbox);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

void interpolate_array_nearest(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target)
{
  if (!validate_non_empty(source) || !validate_non_empty(target)) return;

  glm::vec4 bbox_target_mod = helper_transform_bbox(bbox_source, bbox_target);

  // compute
  auto run = clwrapper::Run("interpolate_array_nearest");

  run.bind_imagef("source", source.vector, source.shape.x, source.shape.y);
  run.bind_imagef("target",
                  target.vector,
                  target.shape.x,
                  target.shape.y,
                  true);

  run.bind_arguments(target.shape.x, target.shape.y, bbox_target_mod);

  run.execute({target.shape.x, target.shape.y});

  run.read_imagef("target");
}

} // namespace hmap::gpu
