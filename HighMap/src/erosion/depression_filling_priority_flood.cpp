/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/drainage_basin_cell_based.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

void depression_filling_priority_flood(Array &z,
                                       bool   apply_post_filter,
                                       bool   outflow_left,
                                       bool   outflow_right,
                                       bool   outflow_bottom,
                                       bool   outflow_top)
{
  if (!validate_non_empty(z)) return;

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
}

} // namespace hmap

namespace hmap::va
{

VirtualArray depression_filling_priority_flood(const VirtualArray &z,
                                               bool          apply_post_filter,
                                               VirtualArray *p_fill_map,
                                               const ComputeMode &cm)
{
  VirtualArray out;
  if (z.empty()) return out;

  out.copy_from(z, cm);

  // iterative multi-pass tile filling
  glm::ivec2 tiling = out.get_max_tiles();
  int        nit = std::max(tiling.x, tiling.y);

  for (int it = 0; it < nit; ++it)
  {
    hmap::for_each_tile(
        {},
        {&out},
        [&](std::vector<const hmap::Array *>,
            std::vector<hmap::Array *> p_arrays_out,
            const hmap::TileRegion &)
        {
          auto [pa_out] = unpack<1>(p_arrays_out);

          hmap::depression_filling_priority_flood(*pa_out, apply_post_filter);
        },
        cm);

    out.sync_overlap_buffers(hmap::SyncOperation::Max);
  }

  // compute fill map
  if (p_fill_map)
  {
    p_fill_map->copy_from(z, cm, /* copy_src_data */ false);

    hmap::for_each_tile(
        {&z, &out},
        {p_fill_map},
        [&](std::vector<const hmap::Array *> p_arrays_in,
            std::vector<hmap::Array *>       p_arrays_out,
            const hmap::TileRegion &)
        {
          auto [pa_in, pa_out] = unpack<2>(p_arrays_in);
          auto [pa_fill_map] = unpack<1>(p_arrays_out);

          *pa_fill_map = *pa_out - *pa_in;
        },
        cm);
  }

  return out;
}

} // namespace hmap::va
