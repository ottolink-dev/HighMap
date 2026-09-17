/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/shortest_path.hpp"

namespace hmap
{

Path dijkstra(const Path  &path,
              const Array &array,
              glm::vec4    bbox,
              float        elevation_ratio,
              float        distance_exponent,
              float        upward_penalization,
              Array       *p_mask_nogo)
{
  if (!validate_non_empty(array)) return path;
  if (!validate_min_size(path, 2, "Path points")) return path;
  if (p_mask_nogo && !validate_same_shape(array, *p_mask_nogo)) return path;

  Path new_path = path;

  size_t ks = new_path.is_closed() ? 0 : 1; // trick to handle closed contours
  for (size_t k = 0; k < new_path.size() - ks; k++)
  {
    size_t knext = (k + 1) % new_path.size();

    glm::ivec2 ij_start = glm::ivec2(
        (int)((new_path[k].x - bbox.x) / (bbox.y - bbox.x) *
              (array.shape.x - 1)),
        (int)((new_path[k].y - bbox.z) / (bbox.w - bbox.z) *
              (array.shape.y - 1)));

    glm::ivec2 ij_end = glm::ivec2(
        (int)((new_path[knext].x - bbox.x) / (bbox.y - bbox.x) *
              (array.shape.x - 1)),
        (int)((new_path[knext].y - bbox.z) / (bbox.w - bbox.z) *
              (array.shape.y - 1)));

    std::vector<int> ip, jp;
    find_path_dijkstra(array,
                       ij_start,
                       ij_end,
                       ip,
                       jp,
                       elevation_ratio,
                       distance_exponent,
                       upward_penalization,
                       p_mask_nogo);

    // backup cuurrent edge informations before adding points to this edge
    Point p1 = new_path[k];
    Point p2 = new_path[knext];

    for (size_t r = 1; r < ip.size() - 1; r++)
    {
      float x = (float)ip[r] / (float)(array.shape.x - 1) * (bbox.y - bbox.x) +
                bbox.x;
      float y = (float)jp[r] / (float)(array.shape.y - 1) * (bbox.w - bbox.z) +
                bbox.z;

      // use distance to start and end points to determine value at the added
      // point (barycentric value)
      float d1 = std::hypot(x - p1.x, y - p1.y);
      float d2 = std::hypot(x - p2.x, y - p2.y);
      float v = (d2 * p1.v + d1 * p2.v) / (d1 + d2);

      Point p = Point(x, y, v);
      new_path.insert(new_path.begin() + k + 1, p);
      ++k;
    }
  }

  return new_path;
}

} // namespace hmap
