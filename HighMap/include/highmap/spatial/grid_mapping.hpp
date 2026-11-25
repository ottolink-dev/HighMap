/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file grid_mapping.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Canonical index <-> UV <-> tile coordinate mapping helpers.
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <algorithm>
#include <cmath>

#include <glm/vec2.hpp>

#include "highmap/spatial/uv_rect.hpp"

namespace hmap
{

/**
 * @brief Grid registration policy for mapping discrete indices to [0, 1].
 */
enum class GridRegistration
{
  CellCenter, ///< Cell i covers [i/N, (i+1)/N], center at (i + 0.5)/N.
  CellCorner, ///< Left corner of cell i at i/N.
  Node        ///< Vertex-centered at i/(N-1), pinning endpoints 0 and 1.
};

/**
 * @brief Maps a discrete index i in [0, n) to a UV coordinate in [0, 1].
 */
inline float index_to_uv(int              i,
                         int              n,
                         GridRegistration reg = GridRegistration::CellCenter)
{
  if (n <= 0) return 0.f;
  switch (reg)
  {
  case GridRegistration::CellCenter:
    return (static_cast<float>(i) + 0.5f) / static_cast<float>(n);
  case GridRegistration::CellCorner:
    return static_cast<float>(i) / static_cast<float>(n);
  case GridRegistration::Node:
    return (n > 1) ? static_cast<float>(i) / static_cast<float>(n - 1) : 0.f;
  }
  return 0.f;
}

/**
 * @brief Maps a UV coordinate in [0, 1] to a continuous index in [0, n - 1].
 */
inline float uv_to_continuous_index(
    float            u,
    int              n,
    GridRegistration reg = GridRegistration::CellCenter)
{
  if (n <= 1) return 0.f;
  float idx = 0.f;
  switch (reg)
  {
  case GridRegistration::CellCenter:
    idx = u * static_cast<float>(n) - 0.5f;
    break;
  case GridRegistration::CellCorner: idx = u * static_cast<float>(n); break;
  case GridRegistration::Node: idx = u * static_cast<float>(n - 1); break;
  }
  return std::clamp(idx, 0.f, static_cast<float>(n - 1));
}

/**
 * @brief Computes the exact UV rectangle of a tile (including halo) within a
 * domain of global_shape.
 */
inline UvRect tile_uv_rect(const glm::ivec2 &global_offset,
                           const glm::ivec2 &tile_shape_with_halo,
                           const glm::ivec2 &global_shape,
                           const UvRect     &domain_uv = UvRect::unit())
{
  const float du = (global_shape.x > 0) ? (domain_uv.u_max - domain_uv.u_min) /
                                              static_cast<float>(global_shape.x)
                                        : 0.f;
  const float dv = (global_shape.y > 0) ? (domain_uv.v_max - domain_uv.v_min) /
                                              static_cast<float>(global_shape.y)
                                        : 0.f;

  return {domain_uv.u_min + static_cast<float>(global_offset.x) * du,
          domain_uv.u_min +
              static_cast<float>(global_offset.x + tile_shape_with_halo.x) * du,
          domain_uv.v_min + static_cast<float>(global_offset.y) * dv,
          domain_uv.v_min +
              static_cast<float>(global_offset.y + tile_shape_with_halo.y) *
                  dv};
}

} // namespace hmap
