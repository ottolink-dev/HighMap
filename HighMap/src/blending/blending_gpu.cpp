/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/range.hpp"
#include "highmap/texture.hpp"

namespace hmap::gpu
{

Array blend_gradients(const Array &array1, const Array &array2, int ir)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();

  Array dn1 = hmap::gradient_norm(array1);
  Array dn2 = hmap::gradient_norm(array2);

  gpu::smooth_cpulse(dn1, ir);
  gpu::smooth_cpulse(dn2, ir);

  Array t = hmap::maximum_smooth(dn1, dn2, 0.1f / (float)array1.shape.x);
  remap(t);

  return lerp(array1, array2, t);
}

Array blend_poisson_bf(const Array &array1,
                       const Array &array2,
                       const int    iterations,
                       const Array *p_mask)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Array();
  if (p_mask && !validate_same_shape(array1, *p_mask)) return Array();

  int nx = array1.shape.x;
  int ny = array1.shape.y;

  // Analytically compute optimal SOR relaxation factor omega:
  // rho_jacobi = 0.5 * (cos(pi / nx) + cos(pi / ny))
  // omega_opt  = 2.0 / (1.0 + sqrt(1.0 - rho_jacobi^2))
  float pi = static_cast<float>(M_PI);
  float rho = 0.5f * (std::cos(pi / static_cast<float>(nx)) +
                      std::cos(pi / static_cast<float>(ny)));
  float omega = 2.f / (1.f + std::sqrt(std::max(0.f, 1.f - rho * rho)));

  Array array1_out = p_mask ? lerp(array1, array2, *p_mask) : array1;
  Array delta2(array1.shape);

  // Precompute Laplacian of array2 (RHS)
  auto run_laplacian = clwrapper::Run("poisson_compute_laplacian");
  run_laplacian.bind_buffer<float>("src", array2.vector);
  run_laplacian.bind_buffer<float>("laplacian", delta2.vector);
  run_laplacian.bind_arguments(nx, ny);
  run_laplacian.write_buffer("src");
  run_laplacian.execute({nx, ny});
  run_laplacian.read_buffer("laplacian");

  // Run Red-Black SOR solver
  auto run = clwrapper::Run("blend_poisson_red_black");
  run.bind_buffer<float>("array1", array1_out.vector);
  run.bind_buffer<float>("delta2", delta2.vector);
  helper_bind_optional_buffer(run, "mask", p_mask);
  run.bind_arguments(nx, ny, 0, omega, p_mask ? 1 : 0);

  run.write_buffer("array1");
  run.write_buffer("delta2");

  for (int it = 0; it < iterations; it++)
  {
    // Pass 0: Red cells ((x + y) % 2 == 0)
    run.set_argument(5, 0);
    run.execute({nx, ny});

    // Pass 1: Black cells ((x + y) % 2 == 1)
    run.set_argument(5, 1);
    run.execute({nx, ny});
  }

  run.read_buffer("array1");

  return array1_out;
}

Texture blend_poisson_bf(const Texture &texture1,
                         const Texture &texture2,
                         const int      iterations,
                         const Array   *p_mask)
{
  if (!validate_non_empty(texture1) || !validate_same_shape(texture1, texture2))
    return Texture();
  if (!validate_channels(texture2, texture1.num_channels())) return Texture();
  if (p_mask && !validate_same_shape(texture1, *p_mask)) return Texture();

  std::vector<Array> blended_channels;
  blended_channels.reserve(texture1.num_channels());

  for (size_t c = 0; c < texture1.channels.size(); ++c)
  {
    blended_channels.push_back(blend_poisson_bf(texture1.channels[c],
                                                texture2.channels[c],
                                                iterations,
                                                p_mask));
  }

  return Texture(blended_channels);
}

Array transfer(const Array &source,
               const Array &target,
               int          ir,
               float        amplitude,
               bool         target_prefiltering)
{
  if (!validate_non_empty(source) || !validate_same_shape(source, target))
    return Array();

  // high-pass spatial filter
  Array w = -source;
  gpu::smooth_cpulse(w, ir);
  w += source;

  if (target_prefiltering)
  {
    Array target_f = target;
    gpu::smooth_cpulse(target_f, ir);
    w = target_f + amplitude * w;
  }
  else
  {
    w = target + amplitude * w;
  }

  return w;
}

} // namespace hmap::gpu
