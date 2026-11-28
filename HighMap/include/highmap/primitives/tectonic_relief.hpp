/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file tectonic_relief.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Procedural plate tectonics terrain generator.
 * @copyright Copyright (c) 2025
 */
#pragma once

#include <vector>

#include <glm/vec2.hpp>

#include "highmap/array.hpp"
#include "highmap/spatial/uv_rect.hpp"

namespace hmap
{

/**
 * @brief Represents a single tectonic plate with position, velocity, and crust
 * type.
 */
struct TectonicPlate
{
  glm::vec2 center = {0.5f, 0.5f};
  glm::vec2 velocity = {0.0f, 0.0f};
  bool      is_continental = true;
};

/**
 * @brief Parameters controlling procedural tectonic relief generation.
 */
struct TectonicReliefParams
{
  int       plate_count = 12;
  uint      seed = 0;
  float     continental_ratio = 0.42f;
  float     roughness = 0.5f;
  float     boundary_width = 0.08f;
  glm::vec2 coast_kw = {6.0f, 6.0f};
  glm::vec2 ridge_kw = {12.0f, 12.0f};
  float     relief_amp = 1.0f;
  float     subduction_offset = 0.8f;
  float     fault_amp = 0.15f;
};

/**
 * @brief Scatter tectonic plates on a jittered grid in the unit domain.
 *
 * Plate 0 is always marked continental; remaining plates are marked continental
 * with probability `continental_ratio`.
 *
 * @param  plate_count       Number of plates to generate.
 * @param  seed              Random seed.
 * @param  continental_ratio Probability that each plate (after plate 0) is
 *                           continental.
 * @return                   Vector of scattered tectonic plates.
 */
std::vector<TectonicPlate> scatter_tectonic_plates(
    int   plate_count,
    uint  seed,
    float continental_ratio = 0.42f);

/**
 * @brief Generate a procedural tectonic relief heightmap from explicit plates.
 *
 * @param  shape    Output array shape.
 * @param  plates   Pre-scattered tectonic plates.
 * @param  params   Tectonic relief parameters.
 * @param  p_dx     Optional x domain warp array.
 * @param  p_dy     Optional y domain warp array.
 * @param  p_ranges Optional output mask for mountain ranges / volcanic arcs.
 * @param  p_faults Optional output mask for shear fault ripples.
 * @param  uv_rect  Parametric domain rectangle for tile-consistent evaluation.
 * @return          Generated tectonic relief heightmap.
 */
Array tectonic_relief(glm::ivec2                       shape,
                      const std::vector<TectonicPlate> &plates,
                      const TectonicReliefParams       &params = TectonicReliefParams(),
                      const Array                      *p_dx = nullptr,
                      const Array                      *p_dy = nullptr,
                      Array                            *p_ranges = nullptr,
                      Array                            *p_faults = nullptr,
                      const UvRect                     &uv_rect = UvRect::unit());

/**
 * @brief Generate a procedural tectonic relief heightmap.
 *
 * Scatters `params.plate_count` plates deterministically from `params.seed`
 * and evaluates plate boundaries, subduction zones, rift valleys, mid-ocean
 * ridges, and shear faults.
 *
 * @param  shape    Output array shape.
 * @param  params   Tectonic relief parameters.
 * @param  p_dx     Optional x domain warp array.
 * @param  p_dy     Optional y domain warp array.
 * @param  p_ranges Optional output mask for mountain ranges / volcanic arcs.
 * @param  p_faults Optional output mask for shear fault ripples.
 * @param  uv_rect  Parametric domain rectangle for tile-consistent evaluation.
 * @return          Generated tectonic relief heightmap.
 */
Array tectonic_relief(glm::ivec2                 shape,
                      const TectonicReliefParams &params = TectonicReliefParams(),
                      const Array               *p_dx = nullptr,
                      const Array               *p_dy = nullptr,
                      Array                     *p_ranges = nullptr,
                      Array                     *p_faults = nullptr,
                      const UvRect              &uv_rect = UvRect::unit());

} // namespace hmap
