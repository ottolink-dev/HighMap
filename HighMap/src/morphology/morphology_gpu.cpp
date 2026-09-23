/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"
#include "highmap/morphology.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

Array border(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return array - gpu::erosion(array, ir, kernel_type);
}

Array closing(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return gpu::erosion(gpu::dilation(array, ir, kernel_type), ir, kernel_type);
}

Array closing_by_reconstruction(const Array &array,
                                int          ir,
                                float        k_smooth_max,
                                MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  Array marker = gpu::dilation(array, ir, kernel_type);
  return gpu::reconstruction_by_erosion(marker,
                                        array,
                                        ir,
                                        k_smooth_max,
                                        kernel_type);
}

Array contour_smoothing(const Array &array, int ir, float transition_ratio)
{
  if (!validate_non_empty(array)) return Array();

  Array edt = distance_transform(is_zero(array)) - distance_transform(array);
  gpu::smooth_cpulse(edt, 2 * ir);

  float width = transition_ratio * ir;
  edt /= width;
  clamp(edt, -1.f, 1.f);
  edt = smoothstep3(0.5f * edt + 0.5f);

  return edt;
}

Array dilation(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return gpu::local_max(array, ir, kernel_type);
}

Array dilation_expand_border_only(const Array &array,
                                  int          ir,
                                  MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();

  const glm::ivec2 &shape = array.shape;
  Array             out = gpu::dilation(array, ir, kernel_type);

  // only keep result in the "background" to leave initial vlaues
  // untouched
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      if (array(i, j) != 0.f) out(i, j) = array(i, j);
    }

  return out;
}

Array erosion(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return gpu::local_min(array, ir, kernel_type);
}

Array morphological_black_hat(const Array &array,
                              int          ir,
                              MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return gpu::closing(array, ir, kernel_type) - array;
}

Array morphological_gradient(const Array &array,
                             int          ir,
                             MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  float vmin = array.min();
  return gpu::dilation(array - vmin, ir, kernel_type) -
         gpu::erosion(array - vmin, ir, kernel_type);
}

Array morphological_laplacian(const Array &array,
                              int          ir,
                              MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  float vmin = array.min();
  return gpu::dilation(array - vmin, ir, kernel_type) +
         gpu::erosion(array - vmin, ir, kernel_type) - 2.f * array;
}

Array morphological_top_hat(const Array &array,
                            int          ir,
                            MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return array - gpu::opening(array, ir, kernel_type);
}

Array opening(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return gpu::dilation(gpu::erosion(array, ir, kernel_type), ir, kernel_type);
}

Array opening_by_reconstruction(const Array &array,
                                int          ir,
                                float        k_smooth_min,
                                MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  Array marker = gpu::erosion(array, ir, kernel_type);
  return gpu::reconstruction_by_dilation(marker,
                                         array,
                                         ir,
                                         k_smooth_min,
                                         kernel_type);
}

Array reconstruction_by_dilation(const Array &marker,
                                 const Array &mask,
                                 int          ir,
                                 float        k_smooth_min,
                                 MinMaxKernel kernel_type)
{
  if (!validate_non_empty(marker) || !validate_same_shape(marker, mask))
    return Array();

  constexpr float tol = 1e-6f;
  Array           current = marker;
  Array           next;

  while (true)
  {
    next = gpu::dilation(current, ir, kernel_type);
    next = hmap::minimum_smooth(next, mask, k_smooth_min);

    float diff = 0.f;
    for (int j = 0; j < current.shape.y; ++j)
      for (int i = 0; i < current.shape.x; ++i)
        diff = std::max(diff, std::abs(next(i, j) - current(i, j)));

    std::swap(current, next);
    if (diff < tol) break;
  }
  return current;
}

Array reconstruction_by_erosion(const Array &marker,
                                const Array &mask,
                                int          ir,
                                float        k_smooth_max,
                                MinMaxKernel kernel_type)
{
  if (!validate_non_empty(marker) || !validate_same_shape(marker, mask))
    return Array();

  constexpr float tol = 1e-6f;
  Array           current = marker;
  Array           next;

  while (true)
  {
    next = gpu::erosion(current, ir, kernel_type);
    next = hmap::maximum_smooth(next, mask, k_smooth_max);

    float diff = 0.f;
    for (int j = 0; j < current.shape.y; ++j)
      for (int i = 0; i < current.shape.x; ++i)
        diff = std::max(diff, std::abs(next(i, j) - current(i, j)));

    std::swap(current, next);
    if (diff < tol) break;
  }
  return current;
}

Array relative_distance_from_skeleton(const Array &array,
                                      const Array &skeleton,
                                      int          ir_search,
                                      int          ir_erosion)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, skeleton))
    return Array();

  const glm::ivec2 &shape = array.shape;

  Array border = array - gpu::erosion(array, ir_erosion);
  Array rdist(shape);

  auto run = clwrapper::Run("relative_distance_from_skeleton");

  run.bind_imagef("array", array.vector, shape.x, shape.y);
  run.bind_imagef("sk", skeleton.vector, shape.x, shape.y);
  run.bind_imagef("border", border.vector, shape.x, shape.y);
  run.bind_imagef("rdist", rdist.vector, shape.x, shape.y, true);
  run.bind_arguments(shape.x, shape.y, ir_search);

  run.execute({shape.x, shape.y});

  run.read_imagef("rdist");

  return rdist;
}

Array relative_distance_from_skeleton(const Array &array,
                                      int          ir_search,
                                      bool         zero_at_borders,
                                      int          ir_erosion)
{
  if (!validate_non_empty(array)) return Array();

  Array sk = gpu::skeleton(array, zero_at_borders);
  return gpu::relative_distance_from_skeleton(array, sk, ir_search, ir_erosion);
}

Array skeleton(const Array &array, bool zero_at_borders)
{
  if (!validate_non_empty(array)) return Array();

  Array sk = generate_buffered_array(array, {1, 1, 1, 1});
  set_borders(sk, 0.f, 1);

  const glm::ivec2 &shape_pad = sk.shape;

  Array prev;
  Array diff;

  auto run = clwrapper::Run("thinning");

  run.bind_imagef("in", sk.vector, shape_pad.x, shape_pad.y);
  run.bind_imagef("out", sk.vector, shape_pad.x, sk.shape.y, true);
  run.bind_arguments(shape_pad.x, shape_pad.y, 0);

  do
  {
    prev = sk;

    run.set_argument(4, 0); // pass 1
    run.write_imagef("in");
    run.execute({shape_pad.x, shape_pad.y});
    run.read_imagef("out");

    run.set_argument(4, 1); // pass 2
    run.write_imagef("in");
    run.execute({shape_pad.x, shape_pad.y});
    run.read_imagef("out");

    diff = sk - prev;

  } while (count_non_zero(diff) > 0);

  // remove padding
  sk = sk.extract_slice({1, sk.shape.x - 1, 1, sk.shape.y - 1});

  // set border to zero
  if (zero_at_borders) zeroed_borders(sk);

  return sk;
}

} // namespace hmap::gpu
