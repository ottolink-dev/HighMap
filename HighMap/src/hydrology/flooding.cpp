/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <future>
#include <limits>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/erosion/depression_filling.hpp"
#include "highmap/features.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

Array flooding_uniform_level(const Array &z, float zref)
{
  if (!validate_non_empty(z)) return Array();

  Array water_depth(z.shape, 0.f);

  water_depth = zref - z;
  clamp_min(water_depth, 0.f);

  return water_depth;
}

Array flooding_from_boundaries(const Array      &z,
                               float             zref,
                               bool              from_east,
                               bool              from_west,
                               bool              from_north,
                               bool              from_south,
                               const TileRegion &region,
                               ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return Array();

  // --- Fast path for single tile / non-distributed execution

  if (comm.size() <= 1)
  {
    Array water_depth = Array(z.shape, 0.f);

    if (from_east)
    {
      int   i = z.shape.x - 1;
      float vmin = std::numeric_limits<float>::max();
      int   jmin = 0;

      for (int j = 0; j < z.shape.y - 1; ++j)
      {
        if (z(i, j) < vmin)
        {
          vmin = z(i, j);
          jmin = j;
        }
      }

      water_depth = flooding_from_point(z, i, jmin, zref - z(i, jmin));
    }

    if (from_west)
    {
      int   i = 0;
      float vmin = std::numeric_limits<float>::max();
      int   jmin = 0;

      for (int j = 0; j < z.shape.y - 1; ++j)
      {
        if (z(i, j) < vmin)
        {
          vmin = z(i, j);
          jmin = j;
        }
      }

      water_depth = maximum(water_depth,
                            flooding_from_point(z, i, jmin, zref - z(i, jmin)));
    }

    if (from_north)
    {
      int   j = z.shape.y - 1;
      float vmin = std::numeric_limits<float>::max();
      int   imin = 0;

      for (int i = 0; i < z.shape.x - 1; ++i)
      {
        if (z(i, j) < vmin)
        {
          vmin = z(i, j);
          imin = i;
        }
      }

      water_depth = maximum(water_depth,
                            flooding_from_point(z, imin, j, zref - z(imin, j)));
    }

    if (from_south)
    {
      int   j = 0;
      float vmin = std::numeric_limits<float>::max();
      int   imin = 0;

      for (int i = 0; i < z.shape.x - 1; ++i)
      {
        if (z(i, j) < vmin)
        {
          vmin = z(i, j);
          imin = i;
        }
      }

      water_depth = maximum(water_depth,
                            flooding_from_point(z, imin, j, zref - z(imin, j)));
    }

    return water_depth;
  }

  // --- Multi-tile distributed execution

  int global_nx = region.global_shape.x;
  int global_ny = region.global_shape.y;

  if (global_nx <= 0 || global_ny <= 0)
  {
    hmap::log::error(
        "flooding_from_boundaries: invalid TileRegion global_shape ({}x{})",
        global_nx,
        global_ny);
    return Array();
  }

  bool is_east_tile = region.is_east_boundary();
  bool is_west_tile = region.is_west_boundary();
  bool is_north_tile = region.is_north_boundary();
  bool is_south_tile = region.is_south_boundary();

  Array water_depth(z.shape, 0.f);

  // 1. East boundary: x = global_nx - 1, y in [0, global_ny - 1)
  if (from_east)
  {
    float local_min = std::numeric_limits<float>::max();

    if (is_east_tile)
    {
      int lx = z.shape.x - 1;

      for (int ly = 0; ly < z.shape.y; ++ly)
      {
        int gy = region.local_to_global(lx, ly).y;

        if (gy >= 0 && gy < global_ny - 1)
        {
          if (z(lx, ly) < local_min)
          {
            local_min = z(lx, ly);
          }
        }
      }
    }

    float global_min = comm.all_reduce(local_min, ReduceOp::Min);

    float local_best_y = std::numeric_limits<float>::max();

    if (is_east_tile)
    {
      int lx = z.shape.x - 1;

      for (int ly = 0; ly < z.shape.y; ++ly)
      {
        int gy = region.local_to_global(lx, ly).y;

        if (gy >= 0 && gy < global_ny - 1)
        {
          if (z(lx, ly) == global_min)
          {
            if (float(gy) < local_best_y)
            {
              local_best_y = float(gy);
            }
          }
        }
      }
    }

    float best_gy_f = comm.all_reduce(local_best_y, ReduceOp::Min);
    int   seed_i = global_nx - 1;
    int   seed_j = int(best_gy_f);

    float depth_min = zref - global_min;

    water_depth = flooding_from_point(z,
                                      seed_i,
                                      seed_j,
                                      depth_min,
                                      region,
                                      comm);
  }

  // 2. West boundary: x = 0, y in [0, global_ny - 1)
  if (from_west)
  {
    float local_min = std::numeric_limits<float>::max();

    if (is_west_tile)
    {
      int lx = 0;

      for (int ly = 0; ly < z.shape.y; ++ly)
      {
        int gy = region.local_to_global(lx, ly).y;

        if (gy >= 0 && gy < global_ny - 1)
        {
          if (z(lx, ly) < local_min)
          {
            local_min = z(lx, ly);
          }
        }
      }
    }

    float global_min = comm.all_reduce(local_min, ReduceOp::Min);

    float local_best_y = std::numeric_limits<float>::max();

    if (is_west_tile)
    {
      int lx = 0;

      for (int ly = 0; ly < z.shape.y; ++ly)
      {
        int gy = region.local_to_global(lx, ly).y;

        if (gy >= 0 && gy < global_ny - 1)
        {
          if (z(lx, ly) == global_min)
          {
            if (float(gy) < local_best_y)
            {
              local_best_y = float(gy);
            }
          }
        }
      }
    }

    float best_gy_f = comm.all_reduce(local_best_y, ReduceOp::Min);
    int   seed_i = 0;
    int   seed_j = int(best_gy_f);

    float depth_min = zref - global_min;

    water_depth = maximum(
        water_depth,
        flooding_from_point(z, seed_i, seed_j, depth_min, region, comm));
  }

  // 3. North boundary: y = global_ny - 1, x in [0, global_nx - 1)
  if (from_north)
  {
    float local_min = std::numeric_limits<float>::max();

    if (is_north_tile)
    {
      int ly = z.shape.y - 1;

      for (int lx = 0; lx < z.shape.x; ++lx)
      {
        int gx = region.local_to_global(lx, ly).x;

        if (gx >= 0 && gx < global_nx - 1)
        {
          if (z(lx, ly) < local_min)
          {
            local_min = z(lx, ly);
          }
        }
      }
    }

    float global_min = comm.all_reduce(local_min, ReduceOp::Min);

    float local_best_x = std::numeric_limits<float>::max();

    if (is_north_tile)
    {
      int ly = z.shape.y - 1;

      for (int lx = 0; lx < z.shape.x; ++lx)
      {
        int gx = region.local_to_global(lx, ly).x;

        if (gx >= 0 && gx < global_nx - 1)
        {
          if (z(lx, ly) == global_min)
          {
            if (float(gx) < local_best_x)
            {
              local_best_x = float(gx);
            }
          }
        }
      }
    }

    float best_gx_f = comm.all_reduce(local_best_x, ReduceOp::Min);
    int   seed_i = int(best_gx_f);
    int   seed_j = global_ny - 1;

    float depth_min = zref - global_min;

    water_depth = maximum(
        water_depth,
        flooding_from_point(z, seed_i, seed_j, depth_min, region, comm));
  }

  // 4. South boundary: y = 0, x in [0, global_nx - 1)
  if (from_south)
  {
    float local_min = std::numeric_limits<float>::max();

    if (is_south_tile)
    {
      int ly = 0;

      for (int lx = 0; lx < z.shape.x; ++lx)
      {
        int gx = region.local_to_global(lx, ly).x;

        if (gx >= 0 && gx < global_nx - 1)
        {
          if (z(lx, ly) < local_min)
          {
            local_min = z(lx, ly);
          }
        }
      }
    }

    float global_min = comm.all_reduce(local_min, ReduceOp::Min);

    float local_best_x = std::numeric_limits<float>::max();

    if (is_south_tile)
    {
      int ly = 0;

      for (int lx = 0; lx < z.shape.x; ++lx)
      {
        int gx = region.local_to_global(lx, ly).x;

        if (gx >= 0 && gx < global_nx - 1)
        {
          if (z(lx, ly) == global_min)
          {
            if (float(gx) < local_best_x)
            {
              local_best_x = float(gx);
            }
          }
        }
      }
    }

    float best_gx_f = comm.all_reduce(local_best_x, ReduceOp::Min);
    int   seed_i = int(best_gx_f);
    int   seed_j = 0;

    float depth_min = zref - global_min;

    water_depth = maximum(
        water_depth,
        flooding_from_point(z, seed_i, seed_j, depth_min, region, comm));
  }

  return water_depth;
}

Array flooding_from_point(const Array      &z,
                          const int         i,
                          const int         j,
                          float             depth_min,
                          const TileRegion &region,
                          ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return Array();

  if (depth_min == std::numeric_limits<float>::max()) depth_min = 0.f;

  std::vector<glm::ivec2> nbrs =
      {{-1, 0}, {0, 1}, {0, -1}, {1, 0}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1}};

  // --- Fast path for single tile / non-distributed execution

  if (comm.size() <= 1)
  {
    Array water_depth(z.shape, 0.f);

    if (i < 0 || i >= z.shape.x || j < 0 || j >= z.shape.y)
    {
      return water_depth;
    }

    float                   z_level = z(i, j) + depth_min;
    std::vector<glm::ivec2> queue = {{i, j}};

    while (!queue.empty())
    {
      glm::ivec2 ij = queue.back();
      queue.pop_back();

      for (const auto &idx : nbrs)
      {
        glm::ivec2 pq = ij + idx;

        if (pq.x >= 0 && pq.x < z.shape.x && pq.y >= 0 && pq.y < z.shape.y)
        {
          float dz = z_level - z(pq.x, pq.y);

          if (z(pq.x, pq.y) < z_level && dz > water_depth(pq.x, pq.y))
          {
            water_depth(pq.x, pq.y) = dz;
            queue.push_back({pq.x, pq.y});
          }
        }
      }
    }

    return water_depth;
  }

  // --- Multi-tile distributed execution

  if (region.tile_shape.x <= 0 || region.tile_shape.y <= 0)
  {
    hmap::log::error(
        "flooding_from_point: invalid TileRegion tile_shape ({}x{})",
        region.tile_shape.x,
        region.tile_shape.y);
    return Array();
  }

  glm::ivec2 local_seed = region.global_to_local(i, j);

  bool owns_seed = (local_seed.x >= 0 && local_seed.x < z.shape.x &&
                    local_seed.y >= 0 && local_seed.y < z.shape.y);

  float local_z_level = owns_seed ? (z(local_seed.x, local_seed.y) + depth_min)
                                  : -std::numeric_limits<float>::infinity();

  float z_level = comm.all_reduce(local_z_level, ReduceOp::Max);

  Array water_depth(z.shape, 0.f);

  if (std::isinf(z_level) || std::isnan(z_level))
  {
    return water_depth;
  }

  std::vector<glm::ivec2> queue;
  if (owns_seed)
  {
    queue.push_back(local_seed);
  }

  int side = std::max(1, int(std::ceil(std::sqrt(comm.size()))));
  int n_sweeps = std::max(1, 2 * side - 1);

  for (int sweep = 0; sweep < n_sweeps; ++sweep)
  {
    while (!queue.empty())
    {
      glm::ivec2 ij = queue.back();
      queue.pop_back();

      for (const auto &idx : nbrs)
      {
        glm::ivec2 pq = ij + idx;

        if (pq.x >= 0 && pq.x < z.shape.x && pq.y >= 0 && pq.y < z.shape.y)
        {
          float dz = z_level - z(pq.x, pq.y);

          if (z(pq.x, pq.y) < z_level && dz > water_depth(pq.x, pq.y))
          {
            water_depth(pq.x, pq.y) = dz;
            queue.push_back({pq.x, pq.y});
          }
        }
      }
    }

    comm.sync_overlap_buffers(water_depth, region, SyncOperation::Max);

    if (sweep < n_sweeps - 1)
    {
      for (int y = 0; y < z.shape.y; ++y)
        for (int x = 0; x < z.shape.x; ++x)
        {
          if (water_depth(x, y) > 0.f)
          {
            for (const auto &idx : nbrs)
            {
              glm::ivec2 pq = glm::ivec2(x, y) + idx;

              if (pq.x >= 0 && pq.x < z.shape.x && pq.y >= 0 &&
                  pq.y < z.shape.y)
              {
                float dz = z_level - z(pq.x, pq.y);

                if (z(pq.x, pq.y) < z_level && dz > water_depth(pq.x, pq.y))
                {
                  queue.push_back({x, y});
                  break;
                }
              }
            }
          }
        }
    }
  }

  // --- Threshold small residuals to avoid precision noise in buffer zones

  constexpr float eps = 1e-4f;
  chop(water_depth, eps);

  return water_depth;
}

Array flooding_from_point(const Array            &z,
                          const std::vector<int> &i,
                          const std::vector<int> &j,
                          float                   depth_min,
                          const TileRegion       &region,
                          ICommunicator          &comm)
{
  if (!validate_non_empty(z)) return Array();
  if (i.size() != j.size())
  {
    hmap::log::warn("Vector size mismatch: i is {}, j is {}",
                    i.size(),
                    j.size());
    return Array();
  }

  Array water_depth(z.shape, 0.f);

  for (size_t k = 0; k < i.size(); k++)
  {
    water_depth = maximum(
        water_depth,
        flooding_from_point(z, i[k], j[k], depth_min, region, comm));
  }

  return water_depth;
}

Array flooding_lake_system(const Array      &z,
                           float             surface_threshold,
                           const TileRegion &region,
                           ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return Array();

  Array water_depth = z;

  // use a depression filling algo to get the lake zones and depths
  depression_filling_priority_flood(water_depth,
                                    false,
                                    true,
                                    true,
                                    true,
                                    true,
                                    region,
                                    comm);

  // compute water depth and threshold negligible residual precision noise
  water_depth -= z;

  constexpr float eps = 1e-4f;
  chop(water_depth, eps);

  // use a connected components analysis to remove small spots if
  // requested
  if (surface_threshold > 0.f)
  {
    Array labels = connected_components(water_depth, surface_threshold);
    for (int j = 0; j < z.shape.y; j++)
      for (int i = 0; i < z.shape.x; i++)
        if (labels(i, j) == 0.f) water_depth(i, j) = 0.f;
  }

  return water_depth;
}

} // namespace hmap
