/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/sdf.hpp"

namespace hmap::gpu
{

Array generate_riverbed(const Path &path,
                        glm::ivec2  shape,
                        glm::vec4   bbox,
                        bool        bezier_smoothing,
                        float       depth_start,
                        float       depth_end,
                        float       slope_start,
                        float       slope_end,
                        float       shape_exponent_start,
                        float       shape_exponent_end,
                        float       k_smoothing,
                        int         post_filter_ir,
                        Array      *p_noise_x,
                        Array      *p_noise_y,
                        Array      *p_noise_r)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_min_size(path, 2, "Path points")) return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array();

  Array sdf;
  if (bezier_smoothing)
    sdf = gpu::sdf_2d_polyline_bezier(path, shape, bbox, p_noise_x, p_noise_y);
  else
    sdf = gpu::sdf_2d_polyline(path, shape, bbox, p_noise_x, p_noise_y);

  Array dz(shape);

  std::vector<float> xp = path.get_x();
  std::vector<float> yp = path.get_y();

  // kernel
  auto run = clwrapper::Run("generate_riverbed");

  run.bind_buffer<float>("sdf", sdf.vector);
  run.bind_buffer<float>("dz", dz.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "noise_r", p_noise_r);
  run.bind_buffer<float>("xp", xp);
  run.bind_buffer<float>("yp", yp);

  run.bind_arguments(shape.x,
                     shape.y,
                     (int)xp.size(),
                     depth_start,
                     depth_end,
                     slope_start,
                     slope_end,
                     shape_exponent_start,
                     shape_exponent_end,
                     k_smoothing,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_noise_r ? 1 : 0,
                     bbox);

  run.write_buffer("sdf");
  run.write_buffer("xp");
  run.write_buffer("yp");

  run.execute({shape.x, shape.y});

  run.read_buffer("dz");

  if (post_filter_ir > 0) gpu::smooth_cpulse(dz, post_filter_ir);

  return dz;
}

} // namespace hmap::gpu
