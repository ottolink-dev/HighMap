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

namespace hmap
{

/**
 * @brief Resample an array to a target array using bicubic interpolation.
 *
 * @param source         Source array.
 * @param target         Target array.
 * @param endpoint       Whether the upper domain endpoint is included.
 * @param pixel_centered Whether coordinates are pixel-centered.
 */
void interpolate_array_bicubic(const Array &source,
                               Array       &target,
                               bool         endpoint = false,
                               bool         pixel_centered = true);

/**
 * @brief Resample a source array into a target array using bicubic
 * interpolation between two local axis-aligned bounding boxes.
 *
 * @param source         Source array.
 * @param target         Target array.
 * @param bbox_source    Source local bounding box {xmin, xmax, ymin, ymax}.
 * @param bbox_target    Target local bounding box {xmin, xmax, ymin, ymax}.
 * @param endpoint       Whether the upper domain endpoint is included.
 * @param pixel_centered Whether coordinates are pixel-centered.
 */
void interpolate_array_bicubic(const Array     &source,
                               Array           &target,
                               const UvRect &uv_rect_source,
                               const UvRect &uv_rect_target,
                               bool             endpoint = false,
                               bool             pixel_centered = true);

/**
 * @brief Resample an array to a target array using bilinear interpolation.
 *
 * @param source         Source array.
 * @param target         Target array.
 * @param endpoint       Whether the upper domain endpoint is included.
 * @param pixel_centered Whether coordinates are pixel-centered.
 */
void interpolate_array_bilinear(const Array &source,
                                Array       &target,
                                bool         endpoint = false,
                                bool         pixel_centered = true);

/**
 * @brief Resample a source array into a target array using bilinear
 * interpolation between two local axis-aligned bounding boxes.
 *
 * @param source         Source array.
 * @param target         Target array.
 * @param bbox_source    Source local bounding box {xmin, xmax, ymin, ymax}.
 * @param bbox_target    Target local bounding box {xmin, xmax, ymin, ymax}.
 * @param endpoint       Whether the upper domain endpoint is included.
 * @param pixel_centered Whether coordinates are pixel-centered.
 */
void interpolate_array_bilinear(const Array     &source,
                                Array           &target,
                                const UvRect &uv_rect_source,
                                const UvRect &uv_rect_target,
                                bool             endpoint = false,
                                bool             pixel_centered = true);

/**
 * @brief Resample an array to a target array using nearest-neighbor
 * interpolation.
 *
 * @param source   Source array.
 * @param target   Target array.
 * @param endpoint Whether the upper domain endpoint is included.
 */
void interpolate_array_nearest(const Array &source,
                               Array       &target,
                               bool         endpoint = false);

/**
 * @brief Resample a source array into a target array using nearest-neighbor
 * interpolation between two local axis-aligned bounding boxes.
 *
 * @param source      Source array.
 * @param target      Target array.
 * @param bbox_source Source local bounding box {xmin, xmax, ymin, ymax}.
 * @param bbox_target Target local bounding box {xmin, xmax, ymin, ymax}.
 * @param endpoint    Whether the upper domain endpoint is included.
 */
void interpolate_array_nearest(const Array     &source,
                               Array           &target,
                               const UvRect &uv_rect_source,
                               const UvRect &uv_rect_target,
                               bool             endpoint = false);

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
                               const UvRect &uv_rect_source,
                               const UvRect &uv_rect_target);

void interpolate_array_bilinear(const Array &source, Array &target);

void interpolate_array_bilinear(const Array     &source,
                                Array           &target,
                                const UvRect &uv_rect_source,
                                const UvRect &uv_rect_target);

void interpolate_array_lagrange(const Array &source, Array &target, int order);

void interpolate_array_nearest(const Array &source, Array &target);

void interpolate_array_nearest(const Array     &source,
                               Array           &target,
                               const UvRect &uv_rect_source,
                               const UvRect &uv_rect_target);

} // namespace hmap::gpu
