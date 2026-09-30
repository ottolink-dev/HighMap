R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
__constant const int   d8_di[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
__constant const int   d8_dj[8] = {0, 1, -1, 0, -1, 1, -1, 1};
__constant const float d8_dist[8] =
    {1.f, 1.f, 1.f, 1.f, 1.41421356f, 1.41421356f, 1.41421356f, 1.41421356f};
__constant const float d8_inv_dist[8] =
    {1.f, 1.f, 1.f, 1.f, 0.70710678f, 0.70710678f, 0.70710678f, 0.70710678f};
__constant const int d8_rev[8] = {3, 2, 1, 0, 7, 6, 5, 4};
)""
