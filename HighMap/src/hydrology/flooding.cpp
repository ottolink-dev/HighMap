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

Array flooding_from_boundaries(const Array &z,
                               float        zref,
                               bool         from_east,
                               bool         from_west,
                               bool         from_north,
                               bool         from_south)
{
  if (!validate_non_empty(z)) return Array();

  Array water_depth = Array(z.shape);

  // find lowest points on the boundaries and starts the flooding there
  // std::vector<int> is, js;

  if (from_east)
  {
    int   i = z.shape.x - 1;
    float vmin = std::numeric_limits<float>::max();
    int   jmin = 0;

    for (int j = 0; j < z.shape.y - 1; ++j)
      if (z(i, j) < vmin)
      {
        vmin = z(i, j);
        jmin = j;
      }

    water_depth = flooding_from_point(z, i, jmin, zref - z(i, jmin));
  }

  if (from_west)
  {
    int   i = 0;
    float vmin = std::numeric_limits<float>::max();
    int   jmin = 0;

    for (int j = 0; j < z.shape.y - 1; ++j)
      if (z(i, j) < vmin)
      {
        vmin = z(i, j);
        jmin = j;
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
      if (z(i, j) < vmin)
      {
        vmin = z(i, j);
        imin = i;
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
      if (z(i, j) < vmin)
      {
        vmin = z(i, j);
        imin = i;
      }

    water_depth = maximum(water_depth,
                          flooding_from_point(z, imin, j, zref - z(imin, j)));
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

  int my_nominal_w = (region.halo.x == 0) ? (region.shape.x - region.halo.y)
                                          : 0;
  int my_nominal_h = (region.halo.z == 0) ? (region.shape.y - region.halo.w)
                                          : 0;

  float nominal_w_f = comm.all_reduce(float(my_nominal_w), ReduceOp::Max);
  float nominal_h_f = comm.all_reduce(float(my_nominal_h), ReduceOp::Max);

  int tile_shape_x = int(nominal_w_f);
  int tile_shape_y = int(nominal_h_f);

  int gx_start = region.key.tx * tile_shape_x - region.halo.x;
  int gy_start = region.key.ty * tile_shape_y - region.halo.z;

  int local_seed_x = i - gx_start;
  int local_seed_y = j - gy_start;

  bool owns_seed = (local_seed_x >= 0 && local_seed_x < z.shape.x &&
                    local_seed_y >= 0 && local_seed_y < z.shape.y);

  float local_z_level = owns_seed ? (z(local_seed_x, local_seed_y) + depth_min)
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
    queue.push_back({local_seed_x, local_seed_y});
  }

  int n_sweeps = std::max(1, int(std::ceil(std::sqrt(comm.size()))));

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
          bool is_boundary = (x == 0 && region.halo.x > 0) ||
                             (x == z.shape.x - 1 && region.halo.y > 0) ||
                             (y == 0 && region.halo.z > 0) ||
                             (y == z.shape.y - 1 && region.halo.w > 0);

          if (is_boundary && water_depth(x, y) > 0.f)
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
