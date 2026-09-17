/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "highmap/geometry/point_sampling.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

int convert_length_to_pixel(float x,
                            int   nx,
                            bool  lim_inf,
                            bool  lim_sup,
                            float scale)
{
  int ir = x / scale * static_cast<float>(nx);
  if (lim_inf) ir = std::max(ir, 1);
  if (lim_sup) ir = std::min(ir, nx - 1);
  return ir;
}

void expand_points_domain(std::vector<float> &x,
                          std::vector<float> &y,
                          std::vector<float> &value,
                          glm::vec4           bbox)
{
  size_t n = x.size();
  x.resize(9 * n);
  y.resize(9 * n);
  value.resize(9 * n);

  float lx = bbox.y - bbox.x;
  float ly = bbox.w - bbox.z;

  int kshift = 0;
  for (int i = -1; i < 2; i++)
    for (int j = -1; j < 2; j++)
      if ((i != 0) or (j != 0))
      {
        kshift++;
        for (size_t k = 0; k < n; k++)
        {
          x[k + kshift * n] = x[k] + (float)i * lx;
          y[k + kshift * n] = y[k] + (float)j * ly;
          value[k + kshift * n] = value[k];
        }
      }
}

void expand_points_at_domain_boundaries(std::vector<float> &x,
                                        std::vector<float> &y,
                                        std::vector<float> &value,
                                        glm::vec4           bbox,
                                        float               boundary_value)
{
  expand_points_domain_corners(x, y, value, bbox, boundary_value);

  int npoints = std::max(0, (int)std::sqrt((float)x.size()));

  const std::vector<float> xb = linspace(bbox.x, bbox.y, npoints, false);
  const std::vector<float> yb = linspace(bbox.z, bbox.w, npoints, false);

  for (int i = 0; i < npoints; i++)
  {
    x.push_back(xb[i]);
    y.push_back(bbox.z);

    x.push_back(xb[i]);
    y.push_back(bbox.w);

    x.push_back(bbox.x);
    y.push_back(yb[i]);

    x.push_back(bbox.y);
    y.push_back(yb[i]);

    value.push_back(boundary_value);
    value.push_back(boundary_value);
    value.push_back(boundary_value);
    value.push_back(boundary_value);
  }
}

void expand_points_domain_corners(std::vector<float> &x,
                                  std::vector<float> &y,
                                  std::vector<float> &value,
                                  glm::vec4           bbox,
                                  float               corner_value)
{
  x.push_back(bbox.x);
  x.push_back(bbox.x);
  x.push_back(bbox.y);
  x.push_back(bbox.y);

  y.push_back(bbox.z);
  y.push_back(bbox.w);
  y.push_back(bbox.z);
  y.push_back(bbox.w);

  for (int i = 0; i < 4; i++)
    value.push_back(corner_value);
}

void grid_xy_vector(std::vector<float> &x,
                    std::vector<float> &y,
                    glm::ivec2          shape,
                    glm::vec4           bbox,
                    bool                endpoint)
{
  if (!validate_shape(shape))
  {
    x.clear();
    y.clear();
    return;
  }

  x = linspace(bbox.x, bbox.y, shape.x, endpoint);
  y = linspace(bbox.z, bbox.w, shape.y, endpoint);
}

void rescale_grid_from_unit_square_to_bbox(std::vector<float> &x,
                                           std::vector<float> &y,
                                           glm::vec4           bbox)
{
  size_t nx = x.size();
  size_t ny = y.size();

  if (nx == 1)
  {
    x[0] = bbox.x;
  }
  else if (nx > 1)
  {
    for (size_t i = 0; i < nx; ++i)
    {
      float t = (float)i / ((float)nx - 1.f);
      x[i] = bbox.x + (bbox.y - bbox.x) * t;
    }
  }

  if (ny == 1)
  {
    y[0] = bbox.z;
  }
  else if (ny > 1)
  {
    for (size_t j = 0; j < ny; ++j)
    {
      float t = (float)j / ((float)ny - 1.f);
      y[j] = bbox.z + (bbox.w - bbox.z) * t;
    }
  }
}

void rescale_points_to_unit_square(std::vector<float> &x,
                                   std::vector<float> &y,
                                   glm::vec4           bbox)
{
  float lx = bbox.y - bbox.x;
  float ly = bbox.w - bbox.z;

  if (std::abs(lx) > 1e-9f)
  {
    for (size_t k = 0; k < x.size(); k++)
      x[k] = (x[k] - bbox.x) / lx;
  }
  else
  {
    for (size_t k = 0; k < x.size(); k++)
      x[k] = 0.f;
  }

  if (std::abs(ly) > 1e-9f)
  {
    for (size_t k = 0; k < y.size(); k++)
      y[k] = (y[k] - bbox.z) / ly;
  }
  else
  {
    for (size_t k = 0; k < y.size(); k++)
      y[k] = 0.f;
  }
}

} // namespace hmap
