/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

Array sdf_2d_polyline(const Path  &path,
                      glm::ivec2   shape,
                      glm::vec4    bbox,
                      const Array *p_noise_x,
                      const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  if (path.size() < 2)
  {
    hmap::log::error("at least 2 points needed in the Path to compute the SDF");
    return Array(shape);
  }

  Array sdf2(shape);

  // input points
  std::vector<float> xp = path.get_x();
  std::vector<float> yp = path.get_y();

  // kernel
  auto run = clwrapper::Run("sdf_2d_polyline");

  run.bind_buffer<float>("sdf2", sdf2.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  run.bind_buffer<float>("xp", xp);
  run.bind_buffer<float>("yp", yp);

  run.bind_arguments(shape.x,
                     shape.y,
                     (int)xp.size(),
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_buffer("xp");
  run.write_buffer("yp");

  run.execute({shape.x, shape.y});

  run.read_buffer("sdf2");

  return sdf2;
}

Array sdf_2d_polyline_bezier(const Path  &path,
                             glm::ivec2   shape,
                             glm::vec4    bbox,
                             const Array *p_noise_x,
                             const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  if (path.size() < 3)
  {
    hmap::log::error("at least 3 points needed in the Path to compute the SDF");
    return Array(shape);
  }

  Array sdf2(shape);

  // get points, repeat every 3 points to have a continuous Bezier
  // curve and add a point if the input size is even, to ensure start
  // and end points will be taken into account
  std::vector<float> x_input = path.get_x();
  std::vector<float> y_input = path.get_y();
  std::vector<float> xp;
  std::vector<float> yp;

  if (x_input.size() % 2 == 0)
  {
    // not just in the middle to avoid Bezier degeneration
    size_t k = x_input.size();
    float  new_x = lerp(x_input[k - 1], x_input[k - 2], 0.1f);
    float  new_y = lerp(y_input[k - 1], y_input[k - 2], 0.1f);

    x_input.insert(x_input.end() - 1, new_x);
    y_input.insert(y_input.end() - 1, new_y);
  }

  for (size_t k = 0; k < x_input.size() - 2; k += 2)
  {
    xp.push_back(x_input[k]);
    xp.push_back(x_input[k + 1]);
    xp.push_back(x_input[k + 2]);

    yp.push_back(y_input[k]);
    yp.push_back(y_input[k + 1]);
    yp.push_back(y_input[k + 2]);
  }

  // kernel
  auto run = clwrapper::Run("sdf_2d_polyline_bezier");

  run.bind_buffer<float>("sdf2", sdf2.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  run.bind_buffer<float>("xp", xp);
  run.bind_buffer<float>("yp", yp);

  run.bind_arguments(shape.x,
                     shape.y,
                     (int)xp.size(),
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_buffer("xp");
  run.write_buffer("yp");

  run.execute({shape.x, shape.y});

  run.read_buffer("sdf2");

  return sdf2;
}

} // namespace hmap::gpu
