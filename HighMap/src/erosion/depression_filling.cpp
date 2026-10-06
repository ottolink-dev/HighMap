/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

void depression_filling(Array &z,
                        int    iterations,
                        float  epsilon,
                        bool   outflow_left,
                        bool   outflow_right,
                        bool   outflow_bottom,
                        bool   outflow_top)
{
  if (!validate_non_empty(z)) return;

  const auto &neighbors = neighborhood::MOORE_8;
  const int   nx = z.shape.x;
  const int   ny = z.shape.y;

  Array z_new = z;
  z_new.set_slice({1, nx - 1, 1, ny - 1}, 1e6f);

  if (!outflow_left)
  {
    for (int j = 0; j < ny; j++)
      z_new(0, j) = 1e6f;
  }
  if (!outflow_right)
  {
    for (int j = 0; j < ny; j++)
      z_new(nx - 1, j) = 1e6f;
  }
  if (!outflow_bottom)
  {
    for (int i = 0; i < nx; i++)
      z_new(i, 0) = 1e6f;
  }
  if (!outflow_top)
  {
    for (int i = 0; i < nx; i++)
      z_new(i, ny - 1) = 1e6f;
  }

  int i_min = outflow_left ? 1 : 0;
  int i_max = outflow_right ? nx - 1 : nx;
  int j_min = outflow_bottom ? 1 : 0;
  int j_max = outflow_top ? ny - 1 : ny;

  for (int it = 0; it < iterations; it++)
  {
    for (int j = j_min; j < j_max; j++)
      for (int i = i_min; i < i_max; i++)
      {
        if (z_new(i, j) > z(i, j))
        {
          for (const auto &nbr : neighbors)
          {
            int p = i + nbr.offset.x;
            int q = j + nbr.offset.y;

            if (p < 0 || p >= nx || q < 0 || q >= ny) continue;

            if (z(i, j) >= z_new(p, q) + epsilon * nbr.distance)
            {
              z_new(i, j) = z(i, j);
              break;
            }

            if (z_new(i, j) > z_new(p, q) + epsilon * nbr.distance)
              z_new(i, j) = z_new(p, q) + epsilon * nbr.distance;
          }
        }
      }
  }

  extrapolate_borders(z_new);
  z = z_new;
}

} // namespace hmap
