/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion/depression_filling.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"

namespace hmap
{

void depression_filling(Array            &z,
                        int               iterations,
                        float             epsilon,
                        bool              outflow_left,
                        bool              outflow_right,
                        bool              outflow_bottom,
                        bool              outflow_top,
                        const TileRegion &region,
                        ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return;

  // --- Single-tile execution

  if (comm.size() <= 1)
  {
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
    return;
  }

  // --- Multi-tile distributed execution with communicator synchronization

  const auto &neighbors = neighborhood::MOORE_8;
  const int   nx = z.shape.x;
  const int   ny = z.shape.y;

  // Domain boundary vs inner tile boundary
  bool tile_outflow_l = region.is_west_boundary() ? outflow_left : true;
  bool tile_outflow_r = region.is_east_boundary() ? outflow_right : true;
  bool tile_outflow_b = region.is_south_boundary() ? outflow_bottom : true;
  bool tile_outflow_t = region.is_north_boundary() ? outflow_top : true;

  int i_min = tile_outflow_l ? 1 : 0;
  int i_max = tile_outflow_r ? nx - 1 : nx;
  int j_min = tile_outflow_b ? 1 : 0;
  int j_max = tile_outflow_t ? ny - 1 : ny;

  Array z_new = z;
  z_new.set_slice({1, nx - 1, 1, ny - 1}, 1e6f);

  if (!tile_outflow_l)
  {
    for (int j = 0; j < ny; j++)
      z_new(0, j) = 1e6f;
  }
  if (!tile_outflow_r)
  {
    for (int j = 0; j < ny; j++)
      z_new(nx - 1, j) = 1e6f;
  }
  if (!tile_outflow_b)
  {
    for (int i = 0; i < nx; i++)
      z_new(i, 0) = 1e6f;
  }
  if (!tile_outflow_t)
  {
    for (int i = 0; i < nx; i++)
      z_new(i, ny - 1) = 1e6f;
  }

  int n_tiles = comm.size();
  int side = std::max(1, int(std::ceil(std::sqrt(n_tiles))));
  int n_sweeps = std::max(1, side);
  int it_per_sweep = std::max(1, iterations / n_sweeps);

  for (int sweep = 0; sweep < n_sweeps; ++sweep)
  {
    for (int it = 0; it < it_per_sweep; it++)
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

    comm.sync_overlap_buffers(z_new, region, SyncOperation::Max);
  }

  z = z_new;
}

} // namespace hmap
