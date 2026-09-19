/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/kd_tree.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"

#include "nanoflann.hpp"

namespace hmap
{

Array interpolate2d_idw(glm::ivec2                shape,
                        const std::vector<float> &x,
                        const std::vector<float> &y,
                        const std::vector<float> &values,
                        const Array              *p_noise_x,
                        const Array              *p_noise_y,
                        glm::vec4                 bbox,
                        float                     distance_exp,
                        float                     radius)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_min_size(x, 2, "Point coordinates x") || x.size() != y.size() ||
      x.size() != values.size())
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  // KD-tree
  KDTree tree(x, y);

  // automatically set the search radius
  if (radius == 0.f)
  {
    size_t    k_nbrs = 32;
    glm::vec2 drange = tree.compute_neighbor_distance_range(k_nbrs);
    radius = 2.5f * drange.y;
  }

  const float radius2 = radius * radius;

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

      // use a radius search to avoid hard neighbor switching
      auto matches = tree.radius_search(xi, yi, radius);

      if (matches.empty())
      {
        out(i, j) = 0.f;
        continue;
      }

      float v = 0.f;
      float sum = 0.f;

      for (const auto &res : matches)
      {
        size_t nk = res.first;
        float  dist = res.second;

        if (dist < 1e-10f)
        {
          v = values[nk];
          break;
        }

        float d = dist / radius2;

        float w = 1.f - d;
        w = smoothstep3(w);

        if (distance_exp == 2.f)
        {
          w /= d;
        }
        else if (distance_exp == 1.f)
        {
          w /= std::sqrt(d);
        }
        else
        {
          w /= std::pow(d, 0.5f * distance_exp);
        }

        v += values[nk] * w;
        sum += w;
      }

      if (sum > 0.f) out(i, j) = v / sum;
    }

  return out;
}

} // namespace hmap
