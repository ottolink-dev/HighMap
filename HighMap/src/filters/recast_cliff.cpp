/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cmath>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/functions.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void recast_cliff(Array &array,
                  float  talus,
                  int    ir,
                  float  amplitude,
                  float  gain,
                  int    iterations,
                  Array *p_cliff_mask)
{
  if (!validate_non_empty(array)) return;

  // compute smoothed base elevation if ir > 0
  Array z_smooth = array;
  if (ir > 0) gpu::smooth_cpulse(z_smooth, ir);

  // compute gradient norm on base
  Array gn = gradient_norm(z_smooth);

  // compute slope excess above talus
  Array dn = gn - talus;
  dn *= static_cast<float>(array.shape.x);
  clamp(dn, 0.f, 1.f);
  if (ir > 0)
  {
    gpu::smooth_cpulse(dn, ir);
    clamp(dn, 0.f, 1.f);
  }

  if (p_cliff_mask) *p_cliff_mask = dn;

  // compute gradient field of original heightmap
  Array gx = gradient_x(array);
  Array gy = gradient_y(array);

  // amplify gradient where slope > talus
  Array scale_field = 1.f + amplitude * (gain - 1.f) * dn;
  gx *= scale_field;
  gy *= scale_field;

  // compute divergence of amplified gradient field (Laplacian RHS)
  Array rho = divergence_from_gradients(gx, gy);

  // solve Poisson equation on GPU with Dirichlet boundary conditions
  int   nx = array.shape.x;
  int   ny = array.shape.y;
  float pi = static_cast<float>(M_PI);
  float rho_sor = 0.5f * (std::cos(pi / static_cast<float>(nx)) +
                          std::cos(pi / static_cast<float>(ny)));
  float omega = 2.f / (1.f + std::sqrt(std::max(0.f, 1.f - rho_sor * rho_sor)));

  Array array_out = array;

  auto run = clwrapper::Run("blend_poisson_red_black");
  run.bind_buffer<float>("array1", array_out.vector);
  run.bind_buffer<float>("delta2", rho.vector);
  gpu::helper_bind_optional_buffer(run, "mask", nullptr);
  run.bind_arguments(nx, ny, 0, omega, 0);

  run.write_buffer("array1");
  run.write_buffer("delta2");

  for (int it = 0; it < iterations; it++)
  {
    // red cells
    run.set_argument(5, 0);
    run.execute({nx, ny});

    // black cells
    run.set_argument(5, 1);
    run.execute({nx, ny});
  }

  run.read_buffer("array1");
  array = std::move(array_out);
}

void recast_cliff(Array       &array,
                  float        talus,
                  int          ir,
                  float        amplitude,
                  const Array *p_mask,
                  float        gain,
                  int          iterations,
                  Array       *p_cliff_mask)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(
      array,
      p_mask,
      [&](Array &a) {
        recast_cliff(a, talus, ir, amplitude, gain, iterations, p_cliff_mask);
      });

  if (p_mask && p_cliff_mask) *p_cliff_mask *= *p_mask;
}

void recast_cliff_directional(Array &array,
                              float  talus,
                              int    ir,
                              float  amplitude,
                              float  angle,
                              float  gain,
                              int    iterations,
                              Array *p_cliff_mask)
{
  if (!validate_non_empty(array)) return;

  Array angle_array(array.shape, angle);
  recast_cliff_directional(array,
                           talus,
                           ir,
                           amplitude,
                           angle_array,
                           nullptr,
                           gain,
                           iterations,
                           p_cliff_mask);
}

void recast_cliff_directional(Array       &array,
                              float        talus,
                              int          ir,
                              float        amplitude,
                              float        angle,
                              const Array *p_mask,
                              float        gain,
                              int          iterations,
                              Array       *p_cliff_mask)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  Array angle_array(array.shape, angle);
  recast_cliff_directional(array,
                           talus,
                           ir,
                           amplitude,
                           angle_array,
                           p_mask,
                           gain,
                           iterations,
                           p_cliff_mask);
}

void recast_cliff_directional(Array       &array,
                              float        talus,
                              int          ir,
                              float        amplitude,
                              const Array &angle,
                              float        gain,
                              int          iterations,
                              Array       *p_cliff_mask)
{
  recast_cliff_directional(array,
                           talus,
                           ir,
                           amplitude,
                           angle,
                           nullptr,
                           gain,
                           iterations,
                           p_cliff_mask);
}

void recast_cliff_directional(Array       &array,
                              float        talus,
                              int          ir,
                              float        amplitude,
                              const Array &angle,
                              const Array *p_mask,
                              float        gain,
                              int          iterations,
                              Array       *p_cliff_mask)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, angle)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  if (p_mask)
  {
    apply_with_mask(array,
                    p_mask,
                    [&](Array &a)
                    {
                      recast_cliff_directional(a,
                                               talus,
                                               ir,
                                               amplitude,
                                               angle,
                                               nullptr,
                                               gain,
                                               iterations,
                                               p_cliff_mask);
                    });

    if (p_cliff_mask) *p_cliff_mask *= *p_mask;
    return;
  }

  // convert angle in degrees to radians
  Array alpha = angle * (static_cast<float>(M_PI) / 180.f);

  // compute smoothed base elevation if ir > 0
  Array z_smooth = array;
  if (ir > 0) gpu::smooth_cpulse(z_smooth, ir);

  // compute gradient norm on base
  Array gn = gradient_norm(z_smooth);

  // compute slope excess above talus
  Array dn = gn - talus;
  dn *= static_cast<float>(array.shape.x);
  clamp(dn, 0.f, 1.f);
  if (ir > 0)
  {
    gpu::smooth_cpulse(dn, ir);
    clamp(dn, 0.f, 1.f);
  }

  // orientation scaling
  Array da = gradient_angle(z_smooth);
  da -= alpha;
  da = cos(da);
  clamp_min(da, 0.f);
  if (ir > 0) gpu::smooth_cpulse(da, ir);

  if (p_cliff_mask) *p_cliff_mask = dn * da;

  // compute gradient field of original heightmap
  Array gx = gradient_x(array);
  Array gy = gradient_y(array);

  // amplify gradient where slope > talus and facing target angle
  Array scale_field = 1.f + amplitude * (gain - 1.f) * (dn * da);
  gx *= scale_field;
  gy *= scale_field;

  // compute divergence of amplified gradient field (Laplacian RHS)
  Array rho = divergence_from_gradients(gx, gy);

  // solve Poisson equation on GPU with Dirichlet boundary conditions
  int   nx = array.shape.x;
  int   ny = array.shape.y;
  float pi = static_cast<float>(M_PI);
  float rho_sor = 0.5f * (std::cos(pi / static_cast<float>(nx)) +
                          std::cos(pi / static_cast<float>(ny)));
  float omega = 2.f / (1.f + std::sqrt(std::max(0.f, 1.f - rho_sor * rho_sor)));

  Array array_out = array;

  auto run = clwrapper::Run("blend_poisson_red_black");
  run.bind_buffer<float>("array1", array_out.vector);
  run.bind_buffer<float>("delta2", rho.vector);
  gpu::helper_bind_optional_buffer(run, "mask", nullptr);
  run.bind_arguments(nx, ny, 0, omega, 0);

  run.write_buffer("array1");
  run.write_buffer("delta2");

  for (int it = 0; it < iterations; it++)
  {
    // red cells
    run.set_argument(5, 0);
    run.execute({nx, ny});

    // black cells
    run.set_argument(5, 1);
    run.execute({nx, ny});
  }

  run.read_buffer("array1");
  array = std::move(array_out);
}

} // namespace hmap
