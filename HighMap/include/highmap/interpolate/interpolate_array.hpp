/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file interpolator_array.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief
 *
 * @copyright Copyright (c) 2025 Otto Link
 */
#pragma once
#include "highmap/array.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

class ComputeMode;  // forward decl.
class VirtualArray; // forward decl.

void interpolate_array_bicubic(const Array &source,
                               Array       &target,
                               bool         endpoint = false,
                               bool         pixel_centered = true);

void interpolate_array_bicubic(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target,
                               bool             endpoint = false,
                               bool             pixel_centered = true);

void interpolate_array_bilinear(const Array &source,
                                Array       &target,
                                bool         endpoint = false,
                                bool         pixel_centered = true);

void interpolate_array_bilinear(const Array     &source,
                                Array           &target,
                                const glm::vec4 &bbox_source,
                                const glm::vec4 &bbox_target,
                                bool             endpoint = false,
                                bool             pixel_centered = true);

void interpolate_array_nearest(const Array &source,
                               Array       &target,
                               bool         endpoint = false);

void interpolate_array_nearest(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target,
                               bool             endpoint = false);

// virtual arrays

void flatten_heightmap(VirtualArray       &h_source1,
                       const VirtualArray &h_source2,
                       const CoordFrame   &t_source1,
                       const CoordFrame   &t_source2,
                       const ComputeMode  &cm);

void flatten_heightmap(const std::vector<const VirtualArray *> &h_sources,
                       VirtualArray                            &h_target,
                       const std::vector<const CoordFrame *>   &t_sources,
                       const CoordFrame                        &t_target,
                       const ComputeMode                       &cm);

void interpolate_heightmap(const VirtualArray &h_source,
                           VirtualArray       &h_target,
                           const CoordFrame   &t_source,
                           const CoordFrame   &t_target,
                           const ComputeMode  &cm);

/**
 * @brief Resample an array to a target shape using cubic B-spline interpolation
 * with optional domain noise warping.
 *
 * @param  src       Source array to resample.
 * @param  dst_shape Target output shape.
 * @param  warp      Domain warping amplitude (0 for no warping).
 * @param  seed      Random seed for domain warping noise.
 * @return           Array     Resampled (and optionally warped) array.
 */
Array resample_bicubic_warp(const Array  &src,
                            glm::ivec2    dst_shape,
                            float         warp = 0.f,
                            std::uint32_t seed = 0);

} // namespace hmap

namespace hmap::gpu
{

void interpolate_array_bicubic(const Array &source, Array &target);

void interpolate_array_bicubic(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target);

void interpolate_array_bilinear(const Array &source, Array &target);

void interpolate_array_bilinear(const Array     &source,
                                Array           &target,
                                const glm::vec4 &bbox_source,
                                const glm::vec4 &bbox_target);

void interpolate_array_lagrange(const Array &source, Array &target, int order);

void interpolate_array_nearest(const Array &source, Array &target);

void interpolate_array_nearest(const Array     &source,
                               Array           &target,
                               const glm::vec4 &bbox_source,
                               const glm::vec4 &bbox_target);

} // namespace hmap::gpu
