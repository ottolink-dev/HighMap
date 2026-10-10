/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <future>
#include <vector>

#include "highmap/coord_frame.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate_array.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <format>

namespace hmap
{

void flatten_heightmap(VirtualArray       &h_source1,
                       const VirtualArray &h_source2,
                       const CoordFrame   &t_source1,
                       const CoordFrame   &t_source2,
                       const ComputeMode  &cm)
{
  // work on a copy because of overlapping buffers
  VirtualArray h_source1_cpy;
  h_source1_cpy.copy_from(h_source1, cm);

  auto lambda =
      [&h_source1, &h_source2, &t_source1, &t_source2](float            &cell,
                                                       int               i,
                                                       int               j,
                                                       const TileRegion &region)
  {
    glm::vec2 pos = region.cell_center(i, j);
    glm::vec2 global_pos = t_source1.map_to_global_coords(pos.x, pos.y);

    if (t_source2.is_point_within(global_pos.x, global_pos.y))
    {
      float v_source1 = t_source1.get_heightmap_value_bilinear(h_source1,
                                                               global_pos.x,
                                                               global_pos.y);
      float v_source2 = t_source2.get_heightmap_value_bilinear(h_source2,
                                                               global_pos.x,
                                                               global_pos.y);

      // transition between the two heightmaps based on the
      // distance to the bounding box
      float r = t_source2.normalized_shape_factor(global_pos.x, global_pos.y);

      cell = lerp(v_source1, v_source2, r);
    }
  };

  for_each_cell(h_source1_cpy, lambda, cm);

  h_source1.copy_from(h_source1_cpy, cm);
}

void flatten_heightmap(const std::vector<const VirtualArray *> &h_sources,
                       VirtualArray                            &h_target,
                       const std::vector<const CoordFrame *>   &t_sources,
                       const CoordFrame                        &t_target,
                       const ComputeMode                       &cm)
{
  if (!validate_non_empty(h_sources, "h_sources") ||
      !validate_non_empty(t_sources, "t_sources"))
    return;

  if (h_sources.size() != t_sources.size())
  {
    hmap::log::warn("Size mismatch between h_sources and t_sources");
    return;
  }

  // interpolation only for the first layer
  hmap::interpolate_heightmap(*h_sources[0],
                              h_target,
                              *t_sources[0],
                              t_target,
                              cm);

  // process in-place the remaining layers by flattening them on the
  // current state of the target layer: target <= target & source
  for (size_t k = 1; k < h_sources.size(); ++k)
    flatten_heightmap(h_target, *h_sources[k], t_target, *t_sources[k], cm);
}

void interpolate_heightmap(const VirtualArray &h_source,
                           VirtualArray       &h_target,
                           const CoordFrame   &t_source,
                           const CoordFrame   &t_target,
                           const ComputeMode  &cm)
{
  auto lambda = [&h_source,
                 &t_source,
                 &t_target](float &cell, int i, int j, const TileRegion &region)
  {
    glm::vec2 pos = region.cell_center(i, j);
    glm::vec2 global_pos = t_target.map_to_global_coords(pos.x, pos.y);

    cell = t_source.get_heightmap_value_bilinear(h_source,
                                                 global_pos.x,
                                                 global_pos.y);
  };

  for_each_cell(h_target, lambda, cm);
}

} // namespace hmap
