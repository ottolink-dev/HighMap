/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <stdexcept>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"

namespace hmap
{

Array interpolate2d(glm::ivec2                shape,
                    const std::vector<float> &x,
                    const std::vector<float> &y,
                    const std::vector<float> &values,
                    InterpolationMethod2D     interpolation_method,
                    const Array              *p_noise_x,
                    const Array              *p_noise_y,
                    const UvRect             &uv_rect)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  switch (interpolation_method)
  {
  case InterpolationMethod2D::ITP2D_DELAUNAY:
  {
    return interpolate2d_delaunay(shape,
                                  x,
                                  y,
                                  values,
                                  p_noise_x,
                                  p_noise_y,
                                  uv_rect);
  }

  case InterpolationMethod2D::ITP2D_NEAREST:
  {
    return interpolate2d_nearest(shape,
                                 x,
                                 y,
                                 values,
                                 p_noise_x,
                                 p_noise_y,
                                 uv_rect);
  }

  case InterpolationMethod2D::ITP2D_IDW:
  {
    return interpolate2d_idw(shape,
                             x,
                             y,
                             values,
                             p_noise_x,
                             p_noise_y,
                             uv_rect);
  }

  case InterpolationMethod2D::ITP2D_GAUSSIAN:
  {
    return interpolate2d_gaussian(shape,
                                  x,
                                  y,
                                  values,
                                  p_noise_x,
                                  p_noise_y,
                                  uv_rect);
  }

  case InterpolationMethod2D::ITP2D_NNI:
  {
    return interpolate2d_nni(shape,
                             x,
                             y,
                             values,
                             p_noise_x,
                             p_noise_y,
                             uv_rect);
  }

  case InterpolationMethod2D::ITP2D_DELAUNAY_GRADIENT:
  {
    return interpolate2d_delaunay_gradient(shape,
                                           x,
                                           y,
                                           values,
                                           p_noise_x,
                                           p_noise_y,
                                           uv_rect);
  }

  default: throw std::runtime_error("unknown 2D interpolation method");
  }
}

} // namespace hmap
