/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/terrain_tri_mesh.hpp"

namespace hmap
{

Array interpolate2d_delaunay_gradient(glm::ivec2                shape,
                                      const std::vector<float> &x,
                                      const std::vector<float> &y,
                                      const std::vector<float> &values,
                                      const Array              *p_noise_x,
                                      const Array              *p_noise_y,
                                      const UvRect             &uv_rect,
                                      float                     fill_value,
                                      float gradient_scaling)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_min_size(x, 3, "Point coordinates x") || x.size() != y.size() ||
      x.size() != values.size())
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  // triangulate data
  TerrainTriMesh mesh(x, y, values);
  mesh.compute_gradients();

  // grid interpolation
  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, uv_rect, false);

  Array out(shape);

  int last_tri = 0;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float px = xg[i] + dx;
      float py = yg[j] + dy;

      out(i, j) = mesh.interpolate_z_linear_gradient({px, py},
                                                     last_tri,
                                                     fill_value,
                                                     gradient_scaling);
    }

  return out;
}

} // namespace hmap
