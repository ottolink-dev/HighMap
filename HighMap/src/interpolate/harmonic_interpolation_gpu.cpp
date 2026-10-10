/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"

namespace hmap::gpu
{

Array harmonic_interpolation(const Array &array,
                             const Array &mask_fixed_values,
                             int          iterations_max,
                             float        tolerance,
                             float        omega)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_same_shape(array, mask_fixed_values)) return Array();

  Array u = array;

  // If omega <= 0, analytically compute the optimal SOR relaxation factor
  // based on grid dimensions:
  // rho_jacobi = 0.5 * (cos(pi / nx) + cos(pi / ny))
  // omega_opt  = 2.0 / (1.0 + sqrt(1.0 - rho_jacobi^2))
  if (omega <= 0.f)
  {
    float pi = static_cast<float>(M_PI);
    float rho = 0.5f * (std::cos(pi / static_cast<float>(u.shape.x)) +
                        std::cos(pi / static_cast<float>(u.shape.y)));
    omega = 2.f / (1.f + std::sqrt(std::max(0.f, 1.f - rho * rho)));
  }

  auto run = clwrapper::Run("harmonic_interpolation_red_black_diff");
  std::vector<float> max_diff_buf(1, 0.f);

  run.bind_buffer<float>("u", u.vector);
  run.bind_buffer<float>("mask_fixed_values", mask_fixed_values.vector);
  run.bind_buffer<float>("max_diff", max_diff_buf);
  run.bind_arguments(u.shape.x, u.shape.y, 0, omega, 0);

  run.write_buffer("u");
  run.write_buffer("mask_fixed_values");

  const int check_interval = 16;

  for (int it = 0; it < iterations_max; ++it)
  {
    int track_diff = (tolerance > 0.f) &&
                             ((it % check_interval == check_interval - 1) ||
                              (it == iterations_max - 1))
                         ? 1
                         : 0;

    if (track_diff)
    {
      max_diff_buf[0] = 0.f;
      run.write_buffer("max_diff");
    }

    // Pass 0: Red cells ((i + j) % 2 == 0)
    run.set_argument(5, 0);
    run.set_argument(7, track_diff);
    run.execute({u.shape.x, u.shape.y});

    // Pass 1: Black cells ((i + j) % 2 == 1)
    run.set_argument(5, 1);
    run.set_argument(7, track_diff);
    run.execute({u.shape.x, u.shape.y});

    if (track_diff)
    {
      run.read_buffer("max_diff");
      if (max_diff_buf[0] < tolerance) break;
    }
  }

  run.read_buffer("u");

  return u;
}

Array harmonic_interpolation(const Array &array,
                             const Array &mask_fixed_values,
                             const Array &dx,
                             const Array &dy,
                             int          iterations_max,
                             float        tolerance,
                             float        omega)
{
  if (!validate_non_empty(array)) return Array();
  if (!validate_same_shape(array, mask_fixed_values)) return Array();
  if (!validate_same_shape(array, dx)) return Array();
  if (!validate_same_shape(array, dy)) return Array();

  Array u = array;

  if (omega <= 0.f)
  {
    float pi = static_cast<float>(M_PI);
    float rho = 0.5f * (std::cos(pi / static_cast<float>(u.shape.x)) +
                        std::cos(pi / static_cast<float>(u.shape.y)));
    omega = 2.f / (1.f + std::sqrt(std::max(0.f, 1.f - rho * rho)));
  }

  auto run = clwrapper::Run("harmonic_interpolation_red_black_diff_aniso");
  std::vector<float> max_diff_buf(1, 0.f);

  run.bind_buffer<float>("u", u.vector);
  run.bind_buffer<float>("mask_fixed_values", mask_fixed_values.vector);
  run.bind_buffer<float>("dx", dx.vector);
  run.bind_buffer<float>("dy", dy.vector);
  run.bind_buffer<float>("max_diff", max_diff_buf);
  run.bind_arguments(u.shape.x, u.shape.y, 0, omega, 0);

  run.write_buffer("u");
  run.write_buffer("mask_fixed_values");
  run.write_buffer("dx");
  run.write_buffer("dy");

  const int check_interval = 16;

  for (int it = 0; it < iterations_max; ++it)
  {
    int track_diff = (tolerance > 0.f) &&
                             ((it % check_interval == check_interval - 1) ||
                              (it == iterations_max - 1))
                         ? 1
                         : 0;

    if (track_diff)
    {
      max_diff_buf[0] = 0.f;
      run.write_buffer("max_diff");
    }

    // Pass 0: Red cells ((i + j) % 2 == 0)
    run.set_argument(7, 0);
    run.set_argument(9, track_diff);
    run.execute({u.shape.x, u.shape.y});

    // Pass 1: Black cells ((i + j) % 2 == 1)
    run.set_argument(7, 1);
    run.set_argument(9, track_diff);
    run.execute({u.shape.x, u.shape.y});

    if (track_diff)
    {
      run.read_buffer("max_diff");
      if (max_diff_buf[0] < tolerance) break;
    }
  }

  run.read_buffer("u");

  return u;
}

} // namespace hmap::gpu
