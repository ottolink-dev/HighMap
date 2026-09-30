/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <limits>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/features.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
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

Array flooding_from_point(const Array &z,
                          const int    i,
                          const int    j,
                          float        depth_min)
{
  if (!validate_non_empty(z)) return Array();

  Array water_depth(z.shape, 0.f);

  std::vector<glm::ivec2> nbrs =
      {{-1, 0}, {0, 1}, {0, -1}, {1, 0}, {-1, -1}, {-1, 1}, {1, -1}, {1, 1}};

  if (depth_min == std::numeric_limits<float>::max()) depth_min = 0.f;
  std::vector<glm::ivec2> queue = {{i, j}};

  // loop around the starting point, anything with elevation lower
  // than the reference elevation is water. If not, the cell is
  // outside the "water" mask
  while (queue.size() > 0)
  {
    glm::ivec2 ij = queue.back();
    queue.pop_back();

    for (auto &idx : nbrs)
    {
      glm::ivec2 pq = ij + idx;

      if (pq.x >= 0 && pq.x < z.shape.x && pq.y >= 0 && pq.y < z.shape.y)
      {
        float dz = z(i, j) + depth_min - z(pq.x, pq.y);
        if (z(pq.x, pq.y) < z(i, j) + depth_min && dz > water_depth(pq.x, pq.y))
        {
          water_depth(pq.x, pq.y) = dz;
          queue.push_back({pq.x, pq.y});
        }
      }
    }
  }

  return water_depth;
}

Array flooding_from_point(const Array            &z,
                          const std::vector<int> &i,
                          const std::vector<int> &j,
                          float                   depth_min)
{
  if (!validate_non_empty(z)) return Array();
  if (i.size() != j.size())
  {
    hmap::log::warn("Vector size mismatch: i is {}, j is {}",
                    i.size(),
                    j.size());
    return Array();
  }

  Array water_depth(z.shape);

  for (size_t k = 0; k < i.size(); k++)
    water_depth = maximum(water_depth,
                          flooding_from_point(z, i[k], j[k], depth_min));

  return water_depth;
}

Array flooding_lake_system(const Array &z, float surface_threshold)
{
  if (!validate_non_empty(z)) return Array();

  Array water_depth = z;

  // use a rough depression filling algo to get the lake zones and depths
  depression_filling_priority_flood(water_depth);

  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
      water_depth(i, j) = std::max(0.f, water_depth(i, j) - z(i, j));

  // use a connected components analysis to remove small spots if
  // requested
  if (surface_threshold)
  {
    Array labels = connected_components(water_depth, surface_threshold);
    for (int j = 0; j < z.shape.y; j++)
      for (int i = 0; i < z.shape.x; i++)
        if (labels(i, j) == 0.f) water_depth(i, j) = 0.f;
  }

  return water_depth;
}

} // namespace hmap

namespace hmap::va
{

VirtualArray flooding_lake_system(const VirtualArray &z,
                                  float               surface_threshold,
                                  const ComputeMode  &cm)
{
  if (z.empty()) return VirtualArray();

  VirtualArray water_depth;

  // use depression filling algo to get lake zones and depths
  depression_filling_priority_flood(z, false, &water_depth, cm);

  // use connected components analysis to remove small spots if requested
  if (surface_threshold > 0.f)
  {
    hmap::for_each_tile(
        {},
        {&water_depth},
        [surface_threshold](std::vector<const hmap::Array *>,
                            std::vector<hmap::Array *> p_arrays_out,
                            const hmap::TileRegion &)
        {
          auto [pa_water] = unpack<1>(p_arrays_out);
          Array labels = connected_components(*pa_water, surface_threshold);
          for (int j = 0; j < pa_water->shape.y; j++)
            for (int i = 0; i < pa_water->shape.x; i++)
              if (labels(i, j) == 0.f) (*pa_water)(i, j) = 0.f;
        },
        cm);
  }

  return water_depth;
}

} // namespace hmap::va
