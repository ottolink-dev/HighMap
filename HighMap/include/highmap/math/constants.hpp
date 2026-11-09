/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file constants.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2023
 */
#pragma once

namespace hmap::cte
{

// pi and derivatives
inline constexpr float pi = 3.14159265358979323846f;
inline constexpr float two_pi = 2.0f * pi;
inline constexpr float half_pi = 0.5f * pi;
inline constexpr float inv_pi = 1.0f / pi;

// roots
inline constexpr float sqrt2 = 1.41421356237309504880f;
inline constexpr float inv_sqrt2 = 0.70710678118654752440f;

// angle conversion helpers
inline constexpr float deg_to_rad = pi / 180.0f;
inline constexpr float rad_to_deg = 180.0f / pi;

} // namespace hmap::cte
