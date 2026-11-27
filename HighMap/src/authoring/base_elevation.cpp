/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

Array base_elevation(glm::ivec2                             shape,
                     const std::vector<std::vector<float>> &values,
                     float                                  width_factor,
                     const Array                           *p_noise_x,
                     const Array                           *p_noise_y,
                     const Array                           *p_stretching,
                     const UvRect                          &uv_rect)
{
  if (!validate_shape(shape)) return Array();

  Array array = Array(shape);

  if (values.empty() || values[0].empty()) return array;
  if (p_noise_x && !validate_same_shape(array, *p_noise_x)) return array;
  if (p_noise_y && !validate_same_shape(array, *p_noise_y)) return array;
  if (p_stretching && !validate_same_shape(array, *p_stretching)) return array;

  // get number of control points per direction
  size_t ni = values.size();
  size_t nj = values[0].size();

  // corresponding coordinates
  float              dxc = 1.f / (float)ni;
  float              dyc = 1.f / (float)nj;
  std::vector<float> xc = linspace(0.5f * dxc, 1.f - 0.5f * dxc, ni);
  std::vector<float> yc = linspace(0.5f * dyc, 1.f - 0.5f * dyc, nj);

  // Gaussian half-widths based on the mesh discretization
  float s2x = 2.f / dxc / dxc / width_factor;
  float s2y = 2.f / dyc / dyc / width_factor;

  for (size_t p = 0; p < ni; p++)
    for (size_t q = 0; q < nj; q++)
    {
      auto lambda = [&xc, &yc, &p, &q, &s2x, &s2y](float x_, float y_, float)
      {
        float r2 = s2x * (x_ - xc[p]) * (x_ - xc[p]) +
                   s2y * (y_ - yc[q]) * (y_ - yc[q]);
        return std::exp(-r2);
      };

      Array array_tmp = Array(array.shape);
      fill_array_using_xy_function(array_tmp,
                                   uv_rect,
                                   nullptr,
                                   p_noise_x,
                                   p_noise_y,
                                   p_stretching,
                                   lambda);
      array += values[p][q] * array_tmp;
    }

  return array;
}

} // namespace hmap
