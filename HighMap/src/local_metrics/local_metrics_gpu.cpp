/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <limits>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"

namespace hmap::gpu
{

Array local_aspect_variance(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array alpha = gradient_angle(array);
  return local_variance(alpha, ir);
}

Array local_max(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();

  switch (kernel_type)
  {
  case MinMaxKernel::DISK: return local_max_disk(array, ir);
  case MinMaxKernel::OCTAGON: return local_max_octagon(array, ir);
  case MinMaxKernel::SQUARE: return local_max_square(array, ir);
  default: throw std::invalid_argument("Unknown MinMaxKernel type");
  }
}

Array local_max_disk(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array out(array.shape);

  auto run = clwrapper::Run("local_max");

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", out.vector, array.shape.x, array.shape.y, true);

  run.bind_arguments(array.shape.x, array.shape.y, ir);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return out;
}

Array local_max_octagon(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  // compute axis-aligned (a) and diagonal (b) decomposition radii
  const int b = static_cast<int>(
      std::round((std::sqrt(2.f) - 1.f) * static_cast<float>(ir)));
  const int a = ir - b;

  Array array_out = array;

  auto run = clwrapper::Run("local_max_octagon");

  run.bind_imagef("in", array_out.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", array_out.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y, a, 0);

  // pass 0: horizontal pass
  run.set_argument(4, a);
  run.set_argument(5, 0);
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  // pass 1: vertical pass
  run.set_argument(4, a);
  run.set_argument(5, 1);
  run.write_imagef("in");
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  if (b > 0)
  {
    // pass 2: main diagonal (+1, +1)
    run.set_argument(4, b);
    run.set_argument(5, 2);
    run.write_imagef("in");
    run.execute({array.shape.x, array.shape.y});
    run.read_imagef("out");

    // pass 3: anti diagonal (+1, -1)
    run.set_argument(4, b);
    run.set_argument(5, 3);
    run.write_imagef("in");
    run.execute({array.shape.x, array.shape.y});
    run.read_imagef("out");
  }

  return array_out;
}

Array local_max_square(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = array;

  auto run = clwrapper::Run("local_max_square");

  run.bind_imagef("in", array_out.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", array_out.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y, ir, 0);

  run.set_argument(5, 0); // row pass
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  run.set_argument(5, 1); // col pass
  run.write_imagef("in");
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  return array_out;
}

Array local_mean(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = array;

  auto run = clwrapper::Run("local_mean");

  run.bind_imagef("in", array_out.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", array_out.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y, ir, 0);

  run.set_argument(5, 0); // row pass
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  run.set_argument(5, 1); // col pass
  run.write_imagef("in");
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  return array_out;
}

Array local_min(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();

  switch (kernel_type)
  {
  case MinMaxKernel::DISK: return local_min_disk(array, ir);
  case MinMaxKernel::OCTAGON: return local_min_octagon(array, ir);
  case MinMaxKernel::SQUARE: return local_min_square(array, ir);
  default: throw std::invalid_argument("Unknown MinMaxKernel type");
  }
}

Array local_min_disk(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array out(array.shape);

  auto run = clwrapper::Run("local_min");

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", out.vector, array.shape.x, array.shape.y, true);

  run.bind_arguments(array.shape.x, array.shape.y, ir);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return out;
}

Array local_min_octagon(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  // compute axis-aligned (a) and diagonal (b) decomposition radii
  const int b = static_cast<int>(
      std::round((std::sqrt(2.f) - 1.f) * static_cast<float>(ir)));
  const int a = ir - b;

  Array array_out = array;

  auto run = clwrapper::Run("local_min_octagon");

  run.bind_imagef("in", array_out.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", array_out.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y, a, 0);

  // pass 0: horizontal pass
  run.set_argument(4, a);
  run.set_argument(5, 0);
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  // pass 1: vertical pass
  run.set_argument(4, a);
  run.set_argument(5, 1);
  run.write_imagef("in");
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  if (b > 0)
  {
    // pass 2: main diagonal (+1, +1)
    run.set_argument(4, b);
    run.set_argument(5, 2);
    run.write_imagef("in");
    run.execute({array.shape.x, array.shape.y});
    run.read_imagef("out");

    // pass 3: anti diagonal (+1, -1)
    run.set_argument(4, b);
    run.set_argument(5, 3);
    run.write_imagef("in");
    run.execute({array.shape.x, array.shape.y});
    run.read_imagef("out");
  }

  return array_out;
}

Array local_min_square(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array array_out = array;

  auto run = clwrapper::Run("local_min_square");

  run.bind_imagef("in", array_out.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", array_out.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y, ir, 0);

  run.set_argument(5, 0); // row pass
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  run.set_argument(5, 1); // col pass
  run.write_imagef("in");
  run.execute({array.shape.x, array.shape.y});
  run.read_imagef("out");

  return array_out;
}

Array local_median_deviation(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array mean = gpu::local_mean(array, ir);
  Array med = gpu::median_pseudo(array, ir); // TODO exact
  return abs(mean - med);
}

Array local_relief(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();
  return gpu::local_max(array, ir, kernel_type) -
         gpu::local_min(array, ir, kernel_type);
}

Array local_skewness(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array out(array.shape);

  auto run = clwrapper::Run("local_skewness");

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", out.vector, array.shape.x, array.shape.y, true);

  run.bind_arguments(array.shape.x, array.shape.y, ir);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return out;
}

Array local_variance(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array out(array.shape);

  auto run = clwrapper::Run("local_variance");

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", out.vector, array.shape.x, array.shape.y, true);

  run.bind_arguments(array.shape.x, array.shape.y, ir);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return out;
}

Array local_z_score(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array out(array.shape);

  auto run = clwrapper::Run("local_z_score");

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", out.vector, array.shape.x, array.shape.y, true);

  run.bind_arguments(array.shape.x, array.shape.y, ir);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return out;
}

Array relative_elevation(const Array &array, int ir, MinMaxKernel kernel_type)
{
  if (!validate_non_empty(array)) return Array();

  Array amin = gpu::local_min(array, ir, kernel_type);
  Array amax = gpu::local_max(array, ir, kernel_type);

  return (array - amin) / (amax - amin + std::numeric_limits<float>::min());
}

Array relative_elevation_square_kernel(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array amin = hmap::local_min(array, ir);
  Array amax = hmap::local_max(array, ir);

  return (array - amin) / (amax - amin + std::numeric_limits<float>::min());
}

Array ruggedness(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array rg(array.shape);

  auto run = clwrapper::Run("ruggedness");

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", rg.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y, ir);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return rg;
}

Array rugosity(const Array &z, int ir, bool convex)
{
  if (!validate_non_empty(z)) return Array();

  Array z_avg(z.shape);
  Array z_std(z.shape);
  Array z_skw(z.shape);
  Array zf = z;
  float tol = 1e-30f;

  // use a kernels only for filtering
  gpu::smooth_cpulse(zf, 2 * ir);
  zf = z - zf;
  z_avg = zf;
  gpu::smooth_cpulse(z_avg, ir);
  z_std = (zf - z_avg) * (zf - z_avg);
  gpu::smooth_cpulse(z_std, ir);
  z_skw = (zf - z_avg) * (zf - z_avg) * (zf - z_avg);

  // last part with dedicated kernel
  auto run = clwrapper::Run("rugosity_post");

  run.bind_buffer("z_skw", z_skw.vector);
  run.bind_buffer("z_std", z_std.vector);
  run.bind_arguments(z.shape.x, z.shape.y, tol, convex ? 1 : 0);

  run.write_buffer("z_skw");
  run.write_buffer("z_std");

  run.execute({z.shape.x, z.shape.y});

  run.read_buffer("z_skw");

  return z_skw;
}

Array topographic_position_index(const Array &array, int ir)
{
  if (!validate_non_empty(array)) return Array();

  Array out(array.shape);

  auto run = clwrapper::Run("topographic_position_index");

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", out.vector, array.shape.x, array.shape.y, true);

  run.bind_arguments(array.shape.x, array.shape.y, ir);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");

  return out;
}

} // namespace hmap::gpu
