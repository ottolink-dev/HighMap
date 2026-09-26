/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"

namespace hmap
{

Array interpolate2d_nni(glm::ivec2                shape,
                        const std::vector<float> &x,
                        const std::vector<float> &y,
                        const std::vector<float> &values,
                        const Array              *p_noise_x,
                        const Array              *p_noise_y,
                        glm::vec4                 bbox)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_min_size(x, 3, "Point coordinates x") || x.size() != y.size() ||
      x.size() != values.size())
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox, /* endpoint */ false);

  std::vector<float> xout, yout;
  xout.reserve(shape.x * shape.y);
  yout.reserve(shape.x * shape.y);

  if (p_noise_x || p_noise_y)
  {
    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
        float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

        xout.push_back(xg[i] + dx);
        yout.push_back(yg[j] + dy);
      }
  }
  else
  {
    for (int j = 0; j < shape.y; j++)
      for (int i = 0; i < shape.x; i++)
      {
        xout.push_back(xg[i]);
        yout.push_back(yg[j]);
      }
  }

  NaturalNeighborInterpolator nn;
  nn.setup_output_points(xout, yout);
  nn.build(x, y);

  std::vector<float> result;
  nn.interpolate(values, result);

  if (result.empty()) return Array(shape);

  Array array_out = Array(shape);

  size_t k = 0;
  for (int j = 0; j < shape.y; j++)
    for (int i = 0; i < shape.x; i++)
    {
      array_out(i, j) = result[k];
      k++;
    }

  extrapolate_borders(array_out);

  return array_out;
}

} // namespace hmap
