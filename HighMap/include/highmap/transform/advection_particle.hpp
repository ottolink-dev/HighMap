/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file advection_particle.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for GPU particle-based advection functions.
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include <vector>

#include "highmap/array.hpp"

namespace hmap
{
class Texture;
}

namespace hmap::gpu
{

/**
 * @brief Performs particle-based advection on a scalar field.
 *
 * This function advects particles through a velocity field (or gradient field)
 * derived from the input array @p z. Each particle carries a value from the @p
 * advected_field and evolves over a number of iterations according to advection
 * parameters such as inertia and persistence. Optionally, advection and masking
 * fields can be applied to control which areas are affected.
 *
 * @param  z                 Input scalar field defining the base flow or
 *                           gradient field used for particle advection.
 * @param  advected_field    The scalar field whose values are transported along
 *                           particle trajectories.
 * @param  iterations        Number of advection iterations to perform per
 *                           particle.
 * @param  nparticles        Number of particles to simulate during advection.
 * @param  seed              Random seed used for initializing particle
 *                           positions and stochastic behavior.
 * @param  reverse           If true, reverses the advection direction (useful
 *                           for backtracing or inverse advection).
 * @param  post_filter       If true, applies a smoothing or filtering operation
 *                           after advection to reduce artifacts.
 * @param  advection_length  Maximum advection step length per iteration
 *                           (controls particle displacement scale).
 * @param  value_persistence Controls how much of a particle’s advected value is
 *                           preserved between iterations (1.0 means full
 *                           preservation, lower values introduce decay).
 * @param  inertia           Controls how much a particle’s previous velocity
 *                           affects its next step (simulates momentum).
 * @param  p_advection_mask  Optional mask array controlling where advection
 *                           occurs (nullptr to disable).
 * @param  p_mask            Optional general mask restricting where output
 *                           values are written (nullptr to disable).
 *
 * @return                   An Array representing the advected version of @p
 *                           advected_field after applying particle advection.
 *
 * **Example**
 * @include ex_advection_particle.cpp
 *
 * **Result**
 * @image html ex_advection_particle.png
 */
Array advection_particle(const Array  &z,
                         const Array  &advected_field,
                         int           iterations,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse = false,
                         bool          post_filter = true,
                         float         post_filter_sigma = 0.125f,
                         float         advection_length = 0.1f,
                         float         value_persistence = 0.99f,
                         float         inertia = 0.f,
                         const Array  *p_advection_mask = nullptr,
                         const Array  *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle(
    const Array                      &z,
    const std::vector<const Array *> &advected_fields,
    int                               iterations,
    int                               nparticles,
    std::uint32_t                     seed,
    bool                              reverse = false,
    bool                              post_filter = true,
    float                             post_filter_sigma = 0.125f,
    float                             advection_length = 0.1f,
    float                             value_persistence = 0.99f,
    float                             inertia = 0.f,
    const Array                      *p_advection_mask = nullptr,
    const Array                      *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle(const Array              &z,
                                      const std::vector<Array> &advected_fields,
                                      int                       iterations,
                                      int                       nparticles,
                                      std::uint32_t             seed,
                                      bool                      reverse = false,
                                      bool         post_filter = true,
                                      float        post_filter_sigma = 0.125f,
                                      float        advection_length = 0.1f,
                                      float        value_persistence = 0.99f,
                                      float        inertia = 0.f,
                                      const Array *p_advection_mask = nullptr,
                                      const Array *p_mask = nullptr);

/// @overload
Array advection_particle(const Array  &z,
                         const Array  &advected_field,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse = false,
                         bool          post_filter = true,
                         float         post_filter_sigma = 0.125f,
                         float         advection_length = 0.1f,
                         float         value_persistence = 0.99f,
                         float         inertia = 0.f,
                         const Array  *p_advection_mask = nullptr,
                         const Array  *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle(
    const Array                      &z,
    const std::vector<const Array *> &advected_fields,
    int                               nparticles,
    std::uint32_t                     seed,
    bool                              reverse = false,
    bool                              post_filter = true,
    float                             post_filter_sigma = 0.125f,
    float                             advection_length = 0.1f,
    float                             value_persistence = 0.99f,
    float                             inertia = 0.f,
    const Array                      *p_advection_mask = nullptr,
    const Array                      *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle(const Array              &z,
                                      const std::vector<Array> &advected_fields,
                                      int                       nparticles,
                                      std::uint32_t             seed,
                                      bool                      reverse = false,
                                      bool         post_filter = true,
                                      float        post_filter_sigma = 0.125f,
                                      float        advection_length = 0.1f,
                                      float        value_persistence = 0.99f,
                                      float        inertia = 0.f,
                                      const Array *p_advection_mask = nullptr,
                                      const Array *p_mask = nullptr);

/// @overload
Array advection_particle(const Array  &dx,
                         const Array  &dy,
                         const Array  &advected_field,
                         int           nparticles,
                         std::uint32_t seed,
                         bool          reverse = false,
                         bool          post_filter = true,
                         float         post_filter_sigma = 0.125f,
                         float         advection_length = 0.1f,
                         float         value_persistence = 0.99f,
                         float         inertia = 0.f,
                         const Array  *p_advection_mask = nullptr,
                         const Array  *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle(
    const Array                      &dx,
    const Array                      &dy,
    const std::vector<const Array *> &advected_fields,
    int                               nparticles,
    std::uint32_t                     seed,
    bool                              reverse = false,
    bool                              post_filter = true,
    float                             post_filter_sigma = 0.125f,
    float                             advection_length = 0.1f,
    float                             value_persistence = 0.99f,
    float                             inertia = 0.f,
    const Array                      *p_advection_mask = nullptr,
    const Array                      *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle(const Array              &dx,
                                      const Array              &dy,
                                      const std::vector<Array> &advected_fields,
                                      int                       nparticles,
                                      std::uint32_t             seed,
                                      bool                      reverse = false,
                                      bool         post_filter = true,
                                      float        post_filter_sigma = 0.125f,
                                      float        advection_length = 0.1f,
                                      float        value_persistence = 0.99f,
                                      float        inertia = 0.f,
                                      const Array *p_advection_mask = nullptr,
                                      const Array *p_mask = nullptr);

/// @overload
Texture advection_particle(const Array   &z,
                           const Texture &texture,
                           int            iterations,
                           int            nparticles,
                           std::uint32_t  seed,
                           bool           reverse = false,
                           bool           post_filter = true,
                           float          post_filter_sigma = 0.125f,
                           float          advection_length = 0.1f,
                           float          value_persistence = 0.99f,
                           float          inertia = 0.f,
                           const Array   *p_advection_mask = nullptr,
                           const Array   *p_mask = nullptr);

/// @overload
Texture advection_particle(const Array   &z,
                           const Texture &texture,
                           int            nparticles,
                           std::uint32_t  seed,
                           bool           reverse = false,
                           bool           post_filter = true,
                           float          post_filter_sigma = 0.125f,
                           float          advection_length = 0.1f,
                           float          value_persistence = 0.99f,
                           float          inertia = 0.f,
                           const Array   *p_advection_mask = nullptr,
                           const Array   *p_mask = nullptr);

/**
 * @brief Performs multiscale particle-based advection on a scalar field.
 *
 * This function advects particles through a gradient field derived from the
 * input heightmap @p z across a geometric hierarchy of resolution levels,
 * applying coarse-to-fine delta upsampling with optional domain warping.
 *
 * @param  z                 Input scalar field defining the base flow.
 * @param  advected_field    The scalar field whose values are transported.
 * @param  seed              Random seed for particle positions and warping.
 * @param  steps_per_level   Number of particle steps / iterations per level.
 * @param  particles_ratio   Ratio of particles relative to grid cells.
 * @param  advection_length  Maximum advection step length per iteration.
 * @param  value_persistence Controls how much of the advected value is
 *                           preserved.
 * @param  inertia           Momentum factor for particle velocity.
 * @param  mix               Blend factor for coarse advection deltas into the
 *                           field.
 * @param  warp              Domain warp amplitude during delta upsampling.
 * @param  reverse           If true, reverses the advection direction.
 * @param  post_filter       If true, applies Laplace smoothing post-filter.
 * @param  post_filter_sigma Laplace smoothing sigma parameter.
 * @param  p_advection_mask  Optional mask controlling where advection occurs.
 * @param  p_mask            Optional general mask restricting output values.
 *
 * @return                   An Array representing the advected field.
 */
Array advection_particle_multiscale(
    const Array            &z,
    const Array            &advected_field,
    std::uint32_t           seed,
    const std::vector<int> &steps_per_level = {4, 2, 1},
    float                   particles_ratio = 0.5f,
    float                   advection_length = 0.1f,
    float                   value_persistence = 0.99f,
    float                   inertia = 0.f,
    float                   mix = 1.f,
    float                   warp = 0.f,
    bool                    reverse = false,
    bool                    post_filter = true,
    float                   post_filter_sigma = 0.125f,
    const Array            *p_advection_mask = nullptr,
    const Array            *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle_multiscale(
    const Array                      &z,
    const std::vector<const Array *> &advected_fields,
    std::uint32_t                     seed,
    const std::vector<int>           &steps_per_level = {4, 2, 1},
    float                             particles_ratio = 0.5f,
    float                             advection_length = 0.1f,
    float                             value_persistence = 0.99f,
    float                             inertia = 0.f,
    float                             mix = 1.f,
    float                             warp = 0.f,
    bool                              reverse = false,
    bool                              post_filter = true,
    float                             post_filter_sigma = 0.125f,
    const Array                      *p_advection_mask = nullptr,
    const Array                      *p_mask = nullptr);

/// @overload
std::vector<Array> advection_particle_multiscale(
    const Array              &z,
    const std::vector<Array> &advected_fields,
    std::uint32_t             seed,
    const std::vector<int>   &steps_per_level = {4, 2, 1},
    float                     particles_ratio = 0.5f,
    float                     advection_length = 0.1f,
    float                     value_persistence = 0.99f,
    float                     inertia = 0.f,
    float                     mix = 1.f,
    float                     warp = 0.f,
    bool                      reverse = false,
    bool                      post_filter = true,
    float                     post_filter_sigma = 0.125f,
    const Array              *p_advection_mask = nullptr,
    const Array              *p_mask = nullptr);

/// @overload
Texture advection_particle_multiscale(
    const Array            &z,
    const Texture          &texture,
    std::uint32_t           seed,
    const std::vector<int> &steps_per_level = {4, 2, 1},
    float                   particles_ratio = 0.5f,
    float                   advection_length = 0.1f,
    float                   value_persistence = 0.99f,
    float                   inertia = 0.f,
    float                   mix = 1.f,
    float                   warp = 0.f,
    bool                    reverse = false,
    bool                    post_filter = true,
    float                   post_filter_sigma = 0.125f,
    const Array            *p_advection_mask = nullptr,
    const Array            *p_mask = nullptr);

} // namespace hmap::gpu
