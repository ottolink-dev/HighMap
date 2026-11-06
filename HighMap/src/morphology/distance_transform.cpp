/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"

float f(int i, float gi)
{
  return (float)(i * i) + gi * gi;
}

int sep(int i, int u, float gi, float gu)
{
  return (int)((u * u - i * i + gu * gu - gi * gi) / (2 * (u - i)));
}

namespace hmap
{

Array distance_transform(const Array      &array,
                         bool              return_squared_distance,
                         const TileRegion &region,
                         ICommunicator    &comm)
{
  if (!validate_non_empty(array)) return Array();

  Array dt = Array(array.shape); // output distance
  Array g = Array(array.shape);
  int   ni = array.shape.x;
  int   nj = array.shape.y;

  // Use a multi-tile-safe upper bound so squared distance does not overflow
  // float32 while remaining larger than any global grid dimension
  float inf = (float)((ni + nj) * std::max(1, comm.size()));

  // --- Phase 1: 1D vertical distance transform along y

  for (int i = 0; i < ni; i++)
  {
    // scan 1
    if (array(i, 0) > 0.f)
      g(i, 0) = 0.f;
    else
      g(i, 0) = inf;

    for (int j = 1; j < nj; j++)
    {
      if (array(i, j) > 0.f)
        g(i, j) = 0.f;
      else
        g(i, j) = 1.f + g(i, j - 1);
    }

    // scan 2
    for (int j = nj - 2; j > -1; j--)
      if (g(i, j + 1) < g(i, j)) g(i, j) = 1.f + g(i, j + 1);
  }

  // Synchronize and propagate vertical distances across tile boundaries
  if (comm.size() > 1)
  {
    const int n_sweeps = std::max(1, comm.size() - 1);

    for (int sweep = 0; sweep < n_sweeps; ++sweep)
    {
      comm.sync_overlap_buffers(g, region, SyncOperation::Min);

      for (int i = 0; i < ni; i++)
      {
        for (int j = 1; j < nj; j++)
          g(i, j) = std::min(g(i, j), 1.f + g(i, j - 1));

        for (int j = nj - 2; j > -1; j--)
          g(i, j) = std::min(g(i, j), 1.f + g(i, j + 1));
      }
    }

    comm.sync_overlap_buffers(g, region, SyncOperation::Min);
  }

  // --- Phase 2: 1D parabolic lower envelope along x

  std::vector<int> s(std::max(ni, nj));
  std::vector<int> t(std::max(ni, nj));

  auto run_row_envelope = [&](const Array &g_in, Array &dt_out)
  {
    for (int j = 0; j < nj; j++)
    {
      int q = 0;
      s[0] = 0;
      t[0] = 0;

      // scan 3
      for (int u = 1; u < ni; u++)
      {
        while ((q >= 0) and
               (f(t[q] - s[q], g_in(s[q], j)) > f(t[q] - u, g_in(u, j))))
          q--;

        if (q < 0)
        {
          q = 0;
          s[0] = u;
        }
        else
        {
          int w = 1 + sep(s[q], u, g_in(s[q], j), g_in(u, j));

          if (w < ni)
          {
            q++;
            s[q] = u;
            t[q] = w;
          }
        }
      }

      // scan 4
      for (int u = ni - 1; u > -1; u--)
      {
        dt_out(u, j) = f(u - s[q], g_in(s[q], j));
        if (u == t[q]) q--;
      }
    }
  };

  if (comm.size() <= 1)
  {
    run_row_envelope(g, dt);
  }
  else
  {
    // Compute signed horizontal offset dx(i, j) = u - x_seed to nearest
    // parabola center in row j so 2D parabola centers propagate accurately
    // across East/West tile boundaries
    Array dx_left(array.shape, 0.f);
    Array best_g_left = g;

    Array best_dx(array.shape, 0.f);
    Array best_g = g;

    auto update_horizontal_envelope = [&]()
    {
      for (int j = 0; j < nj; j++)
      {
        for (int u = 0; u < ni; u++)
        {
          float min_d2 = best_dx(u, j) * best_dx(u, j) +
                         best_g(u, j) * best_g(u, j);
          float b_dx = best_dx(u, j);
          float b_g = best_g(u, j);

          for (int i = 0; i < ni; i++)
          {
            float cand_dx = float(u - i) + best_dx(i, j);
            float cand_g = best_g(i, j);
            float cand_d2 = cand_dx * cand_dx + cand_g * cand_g;

            if (cand_d2 < min_d2)
            {
              min_d2 = cand_d2;
              b_dx = cand_dx;
              b_g = cand_g;
            }
          }

          dx_left(u, j) = b_dx;
          best_g_left(u, j) = b_g;
        }

        for (int u = 0; u < ni; u++)
        {
          best_dx(u, j) = dx_left(u, j);
          best_g(u, j) = best_g_left(u, j);
          dt(u, j) = best_dx(u, j) * best_dx(u, j) +
                     best_g(u, j) * best_g(u, j);
        }
      }
    };

    // Track (best_dx, best_g) from local Meijster envelope
    for (int j = 0; j < nj; j++)
    {
      int q = 0;
      s[0] = 0;
      t[0] = 0;

      for (int u = 1; u < ni; u++)
      {
        while ((q >= 0) and (f(t[q] - s[q], g(s[q], j)) > f(t[q] - u, g(u, j))))
          q--;

        if (q < 0)
        {
          q = 0;
          s[0] = u;
        }
        else
        {
          int w = 1 + sep(s[q], u, g(s[q], j), g(u, j));

          if (w < ni)
          {
            q++;
            s[q] = u;
            t[q] = w;
          }
        }
      }

      for (int u = ni - 1; u > -1; u--)
      {
        best_dx(u, j) = float(u - s[q]);
        best_g(u, j) = g(s[q], j);
        dt(u, j) = f(u - s[q], g(s[q], j));
        if (u == t[q]) q--;
      }
    }

    const int n_sweeps = std::max(1, comm.size() - 1);

    for (int sweep = 0; sweep < n_sweeps; ++sweep)
    {
      // Sync best_dx and best_g from West/East using CopyFirst and CopySecond
      Array dx_west = best_dx;
      Array g_west = best_g;
      comm.sync_overlap_buffers(dx_west, region, SyncOperation::CopyFirst);
      comm.sync_overlap_buffers(g_west, region, SyncOperation::CopyFirst);

      Array dx_east = best_dx;
      Array g_east = best_g;
      comm.sync_overlap_buffers(dx_east, region, SyncOperation::CopySecond);
      comm.sync_overlap_buffers(g_east, region, SyncOperation::CopySecond);

      for (int j = 0; j < nj; j++)
        for (int u = 0; u < ni; u++)
        {
          float d2_curr = best_dx(u, j) * best_dx(u, j) +
                          best_g(u, j) * best_g(u, j);
          float d2_w = dx_west(u, j) * dx_west(u, j) +
                       g_west(u, j) * g_west(u, j);
          float d2_e = dx_east(u, j) * dx_east(u, j) +
                       g_east(u, j) * g_east(u, j);

          if (d2_w < d2_curr)
          {
            d2_curr = d2_w;
            best_dx(u, j) = dx_west(u, j);
            best_g(u, j) = g_west(u, j);
          }

          if (d2_e < d2_curr)
          {
            d2_curr = d2_e;
            best_dx(u, j) = dx_east(u, j);
            best_g(u, j) = g_east(u, j);
          }
        }

      update_horizontal_envelope();
    }

    comm.sync_overlap_buffers(dt, region, SyncOperation::Min);
  }

  if (return_squared_distance)
    return dt;
  else
    return sqrt(dt);
}

Array distance_transform_with_closest(const Array     &array,
                                      Mat<glm::ivec2> &closest,
                                      bool             return_squared_distance)
{
  if (!validate_non_empty(array)) return Array();
  if (closest.shape != array.shape) closest = Mat<glm::ivec2>(array.shape);

  Array dt = Array(array.shape); // output distance
  Array g = Array(array.shape);
  int   ni = array.shape.x;
  int   nj = array.shape.y;
  float inf = (float)(ni + nj);

  // closest indices after row pass
  Mat<int> row_closest_i(array.shape);
  Mat<int> row_closest_j(array.shape);

  // phase 1

  for (int i = 0; i < ni; i++)
  {
    row_closest_i(i, 0) = i;
    row_closest_j(i, 0) = 0;

    // scan 1
    if (array(i, 0) > 0.f)
      g(i, 0) = 0.f;
    else
      g(i, 0) = inf;

    for (int j = 1; j < nj; j++)
    {
      if (array(i, j) > 0.f)
      {
        g(i, j) = 0.f;
        row_closest_i(i, j) = i;
        row_closest_j(i, j) = j;
      }
      else
      {
        g(i, j) = 1.f + g(i, j - 1);
        row_closest_i(i, j) = row_closest_i(i, j - 1);
        row_closest_j(i, j) = row_closest_j(i, j - 1);
      }
    }

    // scan 2
    for (int j = nj - 2; j > -1; j--)
      if (g(i, j + 1) < g(i, j))
      {
        g(i, j) = 1.f + g(i, j + 1);
        row_closest_j(i, j) = row_closest_j(i, j + 1);
      }
  }

  // phase 2
  std::vector<int> s(std::max(ni, nj));
  std::vector<int> t(std::max(ni, nj));

  for (int j = 0; j < nj; j++)
  {
    int q = 0;
    s[0] = 0;
    t[0] = 0;

    // scan 3
    for (int u = 1; u < ni; u++)
    {
      while ((q >= 0) and (f(t[q] - s[q], g(s[q], j)) > f(t[q] - u, g(u, j))))
        q--;

      if (q < 0)
      {
        q = 0;
        s[0] = u;
      }
      else
      {
        int w = 1 + sep(s[q], u, g(s[q], j), g(u, j));

        if (w < ni)
        {
          q++;
          s[q] = u;
          t[q] = w;
        }
      }
    }

    // scan 4
    for (int u = ni - 1; u > -1; u--)
    {
      dt(u, j) = f(u - s[q], g(s[q], j));
      closest(u, j) = {row_closest_i(s[q], j), row_closest_j(s[q], j)};

      if (u == t[q]) q--;
    }
  }

  if (return_squared_distance)
    return dt;
  else
    return sqrt(dt);
}

} // namespace hmap
