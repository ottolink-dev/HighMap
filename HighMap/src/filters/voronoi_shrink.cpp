/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

Array voronoi_shrink(const Array  &array,
                     glm::vec2     kw,
                     float         shrink_factor,
                     float         fill_value,
                     std::uint32_t seed,
                     glm::vec2     jitter,
                     float         angle,
                     const Array  *p_mask,
                     const Array  *p_noise_x,
                     const Array  *p_noise_y,
                     glm::vec4     bbox)
{
  if (!validate_non_empty(array)) return Array();
  if (p_mask && !validate_same_shape(array, *p_mask)) return Array(array.shape);
  if (p_noise_x && !validate_same_shape(array, *p_noise_x))
    return Array(array.shape);
  if (p_noise_y && !validate_same_shape(array, *p_noise_y))
    return Array(array.shape);

  const glm::ivec2 shape = array.shape;
  Array            out(shape);

  auto run = clwrapper::Run("voronoi_shrink");

  run.bind_imagef("in", array.vector, shape.x, shape.y);
  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  helper_bind_optional_buffer(run, "mask", p_mask);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(shape.x,
                     shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     jitter,
                     shrink_factor,
                     fill_value,
                     angle,
                     p_mask ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_imagef("in");
  run.execute({shape.x, shape.y});
  run.read_imagef("out");

  return out;
}

Array voronoi_shrink(const Array  &array,
                     float         kw,
                     float         shrink_factor,
                     float         fill_value,
                     std::uint32_t seed,
                     glm::vec2     jitter,
                     float         angle,
                     const Array  *p_mask,
                     const Array  *p_noise_x,
                     const Array  *p_noise_y,
                     glm::vec4     bbox)
{
  return voronoi_shrink(array,
                        {kw, kw},
                        shrink_factor,
                        fill_value,
                        seed,
                        jitter,
                        angle,
                        p_mask,
                        p_noise_x,
                        p_noise_y,
                        bbox);
}

Array voronoi_shrink(const Array &array,
                     const Cloud &cloud,
                     float        shrink_factor,
                     float        fill_value,
                     const Array *p_mask,
                     const Array *p_noise_x,
                     const Array *p_noise_y,
                     glm::vec4    bbox)
{
  if (!validate_non_empty(array)) return Array();
  if (cloud.empty()) return Array(array.shape);
  if (p_mask && !validate_same_shape(array, *p_mask)) return Array(array.shape);
  if (p_noise_x && !validate_same_shape(array, *p_noise_x))
    return Array(array.shape);
  if (p_noise_y && !validate_same_shape(array, *p_noise_y))
    return Array(array.shape);

  const glm::ivec2 shape = array.shape;
  Array            out(shape);

  std::vector<float> xp = cloud.get_x();
  std::vector<float> yp = cloud.get_y();

  auto run = clwrapper::Run("voronoi_shrink_points");

  run.bind_imagef("in", array.vector, shape.x, shape.y);
  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  helper_bind_optional_buffer(run, "mask", p_mask);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  run.bind_buffer<float>("xp", xp);
  run.bind_buffer<float>("yp", yp);

  run.bind_arguments(shape.x,
                     shape.y,
                     (int)xp.size(),
                     shrink_factor,
                     fill_value,
                     p_mask ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_imagef("in");
  run.write_buffer("xp");
  run.write_buffer("yp");

  run.execute({shape.x, shape.y});

  run.read_imagef("out");

  return out;
}

Array voronoi_shrink(const Array              &array,
                     const std::vector<Point> &points,
                     float                     shrink_factor,
                     float                     fill_value,
                     const Array              *p_mask,
                     const Array              *p_noise_x,
                     const Array              *p_noise_y,
                     glm::vec4                 bbox)
{
  return voronoi_shrink(array,
                        Cloud(points),
                        shrink_factor,
                        fill_value,
                        p_mask,
                        p_noise_x,
                        p_noise_y,
                        bbox);
}

Array voronoi_shrink(const Array                  &array,
                     const std::vector<glm::vec2> &points,
                     float                         shrink_factor,
                     float                         fill_value,
                     const Array                  *p_mask,
                     const Array                  *p_noise_x,
                     const Array                  *p_noise_y,
                     glm::vec4                     bbox)
{
  std::vector<Point> pts;
  pts.reserve(points.size());
  for (const auto &p : points)
  {
    pts.emplace_back(p);
  }

  return voronoi_shrink(array,
                        Cloud(std::move(pts)),
                        shrink_factor,
                        fill_value,
                        p_mask,
                        p_noise_x,
                        p_noise_y,
                        bbox);
}

} // namespace hmap::gpu
