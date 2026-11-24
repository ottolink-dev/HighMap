/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/authoring.hpp"
#include "highmap/geometry/point_sampling.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/range.hpp"
#include "highmap/transform.hpp"
#include "highmap/vectors.hpp"

namespace hmap
{

Array stamping(glm::ivec2                shape,
               const std::vector<float> &xr,
               const std::vector<float> &yr,
               const std::vector<float> &zr,
               Array                     kernel,
               int                       kernel_ir,
               bool                      kernel_scale_radius,
               bool                      kernel_scale_amplitude,
               StampingBlendMethod       blend_method,
               std::uint32_t             seed,
               float                     k_smoothing,
               bool                      kernel_flip,
               bool                      kernel_rotate,
               const UvRect             &uv_rect_array)
{
  if (!validate_shape(shape)) return Array();

  Array array = Array(shape);
  if (!validate_non_empty(kernel)) return array;

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  // normalized node coordinates
  std::vector<float> xrs = xr;
  std::vector<float> yrs = yr;
  rescale_points_to_unit_square(xrs, yrs, uv_rect_array);

  // --- define blending function

  std::function<void(float &, float &)> blend_function;

  switch (blend_method)
  {
    // 'va': value array, 'vk': value kernel

  case StampingBlendMethod::ADD:
    blend_function = [](float &va, float &vk) { va += vk; };
    break;

  case StampingBlendMethod::MAXIMUM:
    blend_function = [](float &va, float &vk) { va = std::max(va, vk); };
    break;

  case StampingBlendMethod::MAXIMUM_SMOOTH:
    blend_function = [&k_smoothing](float &va, float &vk)
    { va = maximum_smooth(va, vk, k_smoothing); };
    break;

  case StampingBlendMethod::MINIMUM:
    blend_function = [](float &va, float &vk) { va = std::min(va, vk); };
    break;

  case StampingBlendMethod::MINIMUM_SMOOTH:
    blend_function = [&k_smoothing](float &va, float &vk)
    { va = minimum_smooth(va, vk, k_smoothing); };
    break;

  case StampingBlendMethod::MULTIPLY:
    blend_function = [](float &va, float &vk) { va *= vk; };
    break;

  case StampingBlendMethod::SUBSTRACT:
    blend_function = [](float &va, float &vk) { va -= vk; };
    break;
  }

  // --- do the stamping

  // scale kernel
  glm::ivec2 scaled_kernel_shape = {2 * kernel_ir + 1, 2 * kernel_ir + 1};

  Array kernel_scaled = kernel.resample_to_shape(scaled_kernel_shape);

  // sort points by value to limit the number of a kernel resizing
  std::vector<size_t> ki = argsort(zr);

  for (size_t k : ki)
  {

    // resize kernel (if requested)
    if (kernel_scale_radius)
    {
      glm::ivec2 new_scaled_kernel_shape = {
          std::max(3, (int)(zr[k] * (2 * kernel_ir + 1))),
          std::max(3, (int)(zr[k] * (2 * kernel_ir + 1)))};

      if (new_scaled_kernel_shape != scaled_kernel_shape)
      {
        scaled_kernel_shape = new_scaled_kernel_shape;
        kernel_scaled = kernel.resample_to_shape(scaled_kernel_shape);
      }
    }

    Array kernel_local = kernel_scaled;

    // transform
    if (kernel_flip)
    {
      if (dis(gen) > 0.5f) flip_ud(kernel_local);
      if (dis(gen) > 0.5f) flip_lr(kernel_local);
      if (dis(gen) > 0.5f) rot90(kernel_local);
      if (dis(gen) > 0.5f) kernel_local = transpose(kernel_local);
    }

    // any rotation angle, time consuming and add some scaling
    // distortions to the input kernel

    if (kernel_rotate)
    {
      bool zoom_in = false;
      bool zero_padding = true;
      rotate(kernel_local, 360.f * dis(gen), zoom_in, zero_padding);
    }

    // amplitude
    if (kernel_scale_amplitude) kernel_local *= zr[k];

    // blend (center kernel on point)
    int i0 = (int)(xrs[k] * (shape.x - 1));
    int j0 = (int)(yrs[k] * (shape.y - 1));

    i0 -= (int)(0.5f * scaled_kernel_shape.x);
    j0 -= (int)(0.5f * scaled_kernel_shape.y);

    for (int j = 0; j < scaled_kernel_shape.y; j++)
      for (int i = 0; i < scaled_kernel_shape.x; i++)
      {
        if (i0 + i >= 0 && i0 + i < shape.x && j0 + j >= 0 && j0 + j < shape.y)
          blend_function(array(i0 + i, j0 + j), kernel_local(i, j));
      }
  }

  return array;
}

} // namespace hmap
