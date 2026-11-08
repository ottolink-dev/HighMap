/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/erosion/depression_filling.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/drainage_basin_cell_based.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

void depression_filling_priority_flood(Array            &z,
                                       bool              apply_post_filter,
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
    Array z_bckp = apply_post_filter ? z : Array();

    auto db = DrainageBasinCellBased(z);

    std::vector<glm::ivec2> outlets;
    outlets.reserve(2 * (z.shape.x + z.shape.y));

    for (int j = 0; j < z.shape.y; ++j)
      for (int i = 0; i < z.shape.x; ++i)
      {
        if ((outflow_left && i == 0) || (outflow_right && i == z.shape.x - 1) ||
            (outflow_bottom && j == 0) || (outflow_top && j == z.shape.y - 1))
        {
          outlets.push_back({i, j});
        }
      }

    db.set_outlets(outlets);

    db.compute_receivers_priority_flood();
    db.update_traversals();

    auto upstream_traversals = db.compute_upstream_traversals();

    for (const auto &indices : upstream_traversals)
      for (size_t k = 0; k < indices.size(); ++k)
      {
        const auto &i = indices[k];
        const auto &r = db.receivers(i);

        float new_z = std::max(z(i), z(r));
        z(i) = new_z;
      }

    if (apply_post_filter)
    {
      Array deposition = z - z_bckp;
      laplace(deposition);
      z = z_bckp + deposition;
    }
    return;
  }

  // --- Multi-tile distributed execution with communicator synchronization

  Array z_bckp = apply_post_filter ? z : Array();

  bool tile_outflow_l = region.is_west_boundary() ? outflow_left : true;
  bool tile_outflow_r = region.is_east_boundary() ? outflow_right : true;
  bool tile_outflow_b = region.is_south_boundary() ? outflow_bottom : true;
  bool tile_outflow_t = region.is_north_boundary() ? outflow_top : true;

  int n_tiles = comm.size();
  int side = std::max(1, int(std::ceil(std::sqrt(n_tiles))));
  int n_sweeps = std::max(1, side);

  for (int sweep = 0; sweep < n_sweeps; ++sweep)
  {
    auto db = DrainageBasinCellBased(z);

    std::vector<glm::ivec2> outlets;
    outlets.reserve(2 * (z.shape.x + z.shape.y));

    for (int j = 0; j < z.shape.y; ++j)
      for (int i = 0; i < z.shape.x; ++i)
      {
        if ((tile_outflow_l && i == 0) ||
            (tile_outflow_r && i == z.shape.x - 1) ||
            (tile_outflow_b && j == 0) ||
            (tile_outflow_t && j == z.shape.y - 1))
        {
          outlets.push_back({i, j});
        }
      }

    db.set_outlets(outlets);

    db.compute_receivers_priority_flood();
    db.update_traversals();

    auto upstream_traversals = db.compute_upstream_traversals();

    for (const auto &indices : upstream_traversals)
      for (size_t k = 0; k < indices.size(); ++k)
      {
        const auto &i = indices[k];
        const auto &r = db.receivers(i);

        float new_z = std::max(z(i), z(r));
        z(i) = new_z;
      }

    comm.sync_overlap_buffers(z, region, SyncOperation::Max);
  }

  if (apply_post_filter)
  {
    Array deposition = z - z_bckp;
    laplace(deposition);
    z = z_bckp + deposition;
    comm.sync_overlap_buffers(z, region, SyncOperation::SmoothBlend);
  }
}

} // namespace hmap
