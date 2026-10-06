/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file grid.hpp
 * @copyright Copyright (c) 2023 Otto Link
 */

#pragma once
#include <array>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"

#include <numbers>

namespace hmap
{

// ==========================================================================
//  Neighborhoods
// ==========================================================================

namespace neighborhood
{

struct Neighbor
{
  glm::ivec2 offset;
  float      distance;
  float      inv_distance;
};

struct NeighborWeighted
{
  glm::ivec2 offset;
  float      distance;
  float      inv_distance;
  float      weight;
};

// 4-neighborhood (Von Neumann: N, E, S, W)
inline constexpr std::array<Neighbor, 4> VON_NEUMANN_4 = {{
    {{0, 1}, 1.0f, 1.0f},
    {{1, 0}, 1.0f, 1.0f},
    {{0, -1}, 1.0f, 1.0f},
    {{-1, 0}, 1.0f, 1.0f},
}};

// Moore neighborhood (orthogonal first, then diagonals)
//  6 2 8
//  1 . 4
//  5 3 7
// clang-format off
inline constexpr std::array<Neighbor, 8> MOORE_8 = {{
	{{-1, 0}, 1.0f, 1.0f},
	{{0, 1}, 1.0f, 1.0f},
	{{0, -1}, 1.0f, 1.0f},
	{{1, 0}, 1.0f, 1.0f},
	{{-1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{-1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
}};
// clang-format on

// D8 Flow direction neighborhood (counter-clockwise / standard D8 code
// ordering):
// 5 6 7
// 4 . 0
// 3 2 1
// clang-format off
inline constexpr std::array<Neighbor, 8> D8_FLOW = {{
	{{1, 0}, 1.0f, 1.0f},
	{{1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{0, -1}, 1.0f, 1.0f},
	{{-1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{-1, 0}, 1.0f, 1.0f},
	{{-1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{0, 1}, 1.0f, 1.0f},
	{{1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
}};
// clang-format on

// D-infinity Flow neighborhood (with Quinn effective contour length weights)
// 5 1 7
// 0 . 3
// 4 2 6
// clang-format off
inline constexpr std::array<NeighborWeighted, 8> DINF_8 = {{
	{{-1, 0}, 1.0f, 1.0f, 0.500f},
	{{0, 1}, 1.0f, 1.0f, 0.500f},
	{{0, -1}, 1.0f, 1.0f, 0.500f},
	{{1, 0}, 1.0f, 1.0f, 0.500f},
	{{-1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>, 0.354f},
	{{-1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>, 0.354f},
	{{1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>, 0.354f},
	{{1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>, 0.354f},
}};
// clang-format on

// Geomorphons 8-neighborhood (ray search pattern)
// 7 0 1
// 6 . 2
// 5 4 3
// clang-format off
inline constexpr std::array<Neighbor, 8> GEOMORPHON_8 = {{
	{{-1, 0}, 1.0f, 1.0f},
	{{-1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{0, 1}, 1.0f, 1.0f},
	{{1, 1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{1, 0}, 1.0f, 1.0f},
	{{1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
	{{0, -1}, 1.0f, 1.0f},
	{{-1, -1}, std::numbers::sqrt2_v<float>, 1.0f / std::numbers::sqrt2_v<float>},
}};
// clang-format on

} // namespace neighborhood

// ==========================================================================
//  Functions
// ==========================================================================

/**
 * @brief Converts a length value to a pixel index in a discretized space.
 *
 * This function maps a floating-point length `x` to an integer pixel index,
 * considering the total number of pixels `nx` and a scaling factor `scale`.
 * Optionally, it can enforce lower and upper limits on the output index.
 *
 * @param  x       The length value to be converted.
 * @param  nx      The number of pixels in the discretized space.
 * @param  lim_inf If nonzero, enforces a minimum index of 1.
 * @param  lim_sup If nonzero, enforces a maximum index of `nx - 1`.
 * @param  scale   The scaling factor relating the length to the pixel space.
 * @return         The computed pixel index.
 */
int convert_length_to_pixel(float x,
                            int   nx,
                            bool  lim_inf = true,
                            bool  lim_sup = false,
                            float scale = 1.f);

/**
 * @brief Return x and y coordinates of a regular grid, as two 1D vectors.
 * @param x[out]   Vector x.
 * @param y[out]   Vector y.
 * @param shape    Shape.
 * @param bbox     Bounding box.
 * @param endpoint Include or not the endpoint.
 */
void grid_xy_vector(std::vector<float> &x,
                    std::vector<float> &y,
                    glm::ivec2          shape,
                    glm::vec4           bbox = {0.f, 1.f, 0.f, 1.f},
                    bool                endpoint = false);

/**
 * @brief Rescale 1D grid vectors from unit square [0, 1] to target bounding
 * box.
 * @param x[in, out] Vector x.
 * @param y[in, out] Vector y.
 * @param bbox       Target bounding box.
 */
void rescale_grid_from_unit_square_to_bbox(std::vector<float> &x,
                                           std::vector<float> &y,
                                           glm::vec4           bbox);

} // namespace hmap
