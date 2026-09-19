/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/kd_tree.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

Array interpolate2d_nearest(glm::ivec2                shape,
                            const std::vector<float> &x,
                            const std::vector<float> &y,
                            const std::vector<float> &values,
                            const Array              *p_noise_x,
                            const Array              *p_noise_y,
                            glm::vec4                 bbox)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_min_size(x, 2, "Point coordinates x") || x.size() != y.size() ||
      x.size() != values.size())
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  // KD-tree
  KDTree tree(x, y);

  // interpolation base grid
  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox, /* endpoint */ false);

  // interpolate
  Array out(shape);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;
      float xi = xg[i] + dx;
      float yi = yg[j] + dy;

      out(i, j) = values[tree.nearest(xi, yi)];
    }

  return out;
}

} // namespace hmap
