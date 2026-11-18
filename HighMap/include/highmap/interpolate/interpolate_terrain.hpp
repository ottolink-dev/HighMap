/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file interpolate_terrain.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief World-coordinate (CoordFrame) terrain interpolation and flattening.
 *
 * @copyright Copyright (c) 2025 Otto Link
 */
#pragma once
#include <vector>

#include "highmap/array.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

class ComputeMode; // forward decl.

// --- Array -> Array (world space via CoordFrame)

/**
 * @brief Blend a second terrain array in-place onto a base terrain array using
 * their world-space coordinate frames.
 *
 * @param h_source1 Base terrain array (modified in-place).
 * @param h_source2 Overlay terrain array.
 * @param t_source1 World coordinate frame of @p h_source1.
 * @param t_source2 World coordinate frame of @p h_source2.
 * @param bilinear  Use bilinear interpolation if true, nearest-neighbor
 *                  otherwise.
 */
void flatten_terrain(Array            &h_source1,
                     const Array      &h_source2,
                     const CoordFrame &t_source1,
                     const CoordFrame &t_source2,
                     bool              bilinear = true);

/**
 * @brief Flatten multiple terrain arrays into a target array using their
 * world-space coordinate frames.
 *
 * @param h_sources Ordered list of source terrain arrays.
 * @param h_target  Target terrain array.
 * @param t_sources World coordinate frames corresponding to @p h_sources.
 * @param t_target  World coordinate frame of @p h_target.
 * @param bilinear  Use bilinear interpolation if true, nearest-neighbor
 *                  otherwise.
 */
void flatten_terrain(const std::vector<const Array *>      &h_sources,
                     Array                                 &h_target,
                     const std::vector<const CoordFrame *> &t_sources,
                     const CoordFrame                      &t_target,
                     bool                                   bilinear = true);

/**
 * @brief Project and interpolate a source terrain array onto a target terrain
 * array using their world-space coordinate frames.
 *
 * @param h_source   Source terrain array.
 * @param h_target   Target terrain array.
 * @param t_source   World coordinate frame of @p h_source.
 * @param t_target   World coordinate frame of @p h_target.
 * @param bilinear   Use bilinear interpolation if true, nearest-neighbor
 *                   otherwise.
 * @param fill_value Value assigned to target cells outside @p t_source.
 */
void interpolate_terrain(const Array      &h_source,
                         Array            &h_target,
                         const CoordFrame &t_source,
                         const CoordFrame &t_target,
                         bool              bilinear = true,
                         float             fill_value = 0.f);

// --- VirtualArray -> VirtualArray (world space via CoordFrame)

/**
 * @brief Blend a second terrain VirtualArray in-place onto a base terrain
 * VirtualArray using their world-space coordinate frames.
 *
 * @param h_source1 Base terrain VirtualArray (modified in-place).
 * @param h_source2 Overlay terrain VirtualArray.
 * @param t_source1 World coordinate frame of @p h_source1.
 * @param t_source2 World coordinate frame of @p h_source2.
 * @param cm        Compute mode controlling tiling and CPU distribution.
 * @param bilinear  Use bilinear interpolation if true, nearest-neighbor
 *                  otherwise.
 */
void flatten_terrain(VirtualArray       &h_source1,
                     const VirtualArray &h_source2,
                     const CoordFrame   &t_source1,
                     const CoordFrame   &t_source2,
                     const ComputeMode  &cm,
                     bool                bilinear = true);

/**
 * @brief Flatten multiple terrain VirtualArrays into a target VirtualArray
 * using their world-space coordinate frames.
 *
 * @param h_sources Ordered list of source terrain VirtualArrays.
 * @param h_target  Target terrain VirtualArray.
 * @param t_sources World coordinate frames corresponding to @p h_sources.
 * @param t_target  World coordinate frame of @p h_target.
 * @param cm        Compute mode controlling tiling and CPU distribution.
 * @param bilinear  Use bilinear interpolation if true, nearest-neighbor
 *                  otherwise.
 */
void flatten_terrain(const std::vector<const VirtualArray *> &h_sources,
                     VirtualArray                            &h_target,
                     const std::vector<const CoordFrame *>   &t_sources,
                     const CoordFrame                        &t_target,
                     const ComputeMode                       &cm,
                     bool                                     bilinear = true);

/**
 * @brief Project and interpolate a source terrain VirtualArray onto a target
 * terrain VirtualArray using their world-space coordinate frames.
 *
 * @param h_source   Source terrain VirtualArray.
 * @param h_target   Target terrain VirtualArray.
 * @param t_source   World coordinate frame of @p h_source.
 * @param t_target   World coordinate frame of @p h_target.
 * @param cm         Compute mode controlling tiling and CPU distribution.
 * @param bilinear   Use bilinear interpolation if true, nearest-neighbor
 *                   otherwise.
 * @param fill_value Value assigned to target cells outside @p t_source.
 */
void interpolate_terrain(const VirtualArray &h_source,
                         VirtualArray       &h_target,
                         const CoordFrame   &t_source,
                         const CoordFrame   &t_target,
                         const ComputeMode  &cm,
                         bool                bilinear = true,
                         float               fill_value = 0.f);

} // namespace hmap
