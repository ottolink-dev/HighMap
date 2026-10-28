/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file sdf.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for signed distance functions.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"

namespace hmap
{

/**
 * @brief Computes the signed distance field (SDF) of a 2D polyline.
 *
 * This function calculates the distance of each grid cell in the specified 2D
 * array to the closest segment of a polyline. The polyline is represented by a
 * `Path` object, and optional noise fields can be added to perturb the array
 * coordinates.
 *
 * @param  path      The polyline path represented as a `Path` object containing
 *                   the nodes.
 * @param  shape     The dimensions (width, height) of the output array grid.
 * @param  bbox      The bounding box of the array grid in the format (xmin,
 *                   xmax, ymin, ymax).
 * @param  p_noise_x Pointer to an optional noise array for perturbing
 *                   x-coordinates (can be nullptr).
 * @param  p_noise_y Pointer to an optional noise array for perturbing
 *                   y-coordinates (can be nullptr).
 *
 * @return           Array The resulting signed distance field (SDF) as a 2D
 *                   array.
 *
 * **Example**
 * @include ex_sdf_polyline.cpp
 *
 * **Result**
 * @image html ex_sdf_polyline.png
 */
Array sdf_2d_polyline(const Path  &path,
                      glm::ivec2   shape,
                      glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                      const Array *p_noise_x = nullptr,
                      const Array *p_noise_y = nullptr,
                      bool         use_bvh = true);

Array sdf_2d_polyline_brute_force(const Path  &path,
                                  glm::ivec2   shape,
                                  glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                                  const Array *p_noise_x = nullptr,
                                  const Array *p_noise_y = nullptr);

Array sdf_2d_polyline_bvh(const Path  &path,
                          glm::ivec2   shape,
                          glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                          const Array *p_noise_x = nullptr,
                          const Array *p_noise_y = nullptr);

/*! @brief See hmap::sdf_2d_polyline, with a Bezier approximation of the path.
 */
Array sdf_2d_polyline_bezier(const Path  &path,
                             glm::ivec2   shape,
                             glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                             const Array *p_noise_x = nullptr,
                             const Array *p_noise_y = nullptr,
                             bool         use_bvh = true);

Array sdf_2d_polyline_bezier_brute_force(const Path &path,
                                         glm::ivec2  shape,
                                         glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
                                         const Array *p_noise_x = nullptr,
                                         const Array *p_noise_y = nullptr);

Array sdf_2d_polyline_bezier_bvh(const Path  &path,
                                 glm::ivec2   shape,
                                 glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                                 const Array *p_noise_x = nullptr,
                                 const Array *p_noise_y = nullptr);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::sdf_2d_polyline */
Array sdf_2d_polyline(const Path  &path,
                      glm::ivec2   shape,
                      glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                      const Array *p_noise_x = nullptr,
                      const Array *p_noise_y = nullptr);

/*! @brief See hmap::sdf_2d_polyline_bezier */
Array sdf_2d_polyline_bezier(const Path  &path,
                             glm::ivec2   shape,
                             glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                             const Array *p_noise_x = nullptr,
                             const Array *p_noise_y = nullptr);

} // namespace hmap::gpu