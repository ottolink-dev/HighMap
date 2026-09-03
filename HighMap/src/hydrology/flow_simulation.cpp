/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <cmath>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap::gpu
{

float helper_vflow_compute_adaptive_dt(const Array &h,
                                       float        viscosity,
                                       float        power)
{
  float hmax = h.max();
  float mobility = std::pow(hmax, power) / viscosity;
  float c = 0.2f; // safety factor
  return c / fmax(mobility, 1e-6f);
}

Array flow_simulation(const Array &z,
                      float        water_height,
                      const Array &depth_map,
                      int          iterations,
                      float        dt,
                      bool         flux_diffusion,
                      float        flux_diffusion_strength,
                      float        dry_out_ratio,
                      const Array *p_rain_map,
                      float        rain_rate,
                      float        evap_rate,
                      bool         outflow_boundaries,
                      Array       *p_vel_u,
                      Array       *p_vel_v)
{
  if (!validate_non_empty(z)) return Array();
  if (!validate_same_shape(z, depth_map)) return Array();
  if (p_rain_map && !validate_same_shape(z, *p_rain_map)) return Array();

  const glm::ivec2 shape = z.shape;

  Array d = water_height * depth_map;
  Array fl(shape); // left flux
  Array fr(shape); // right
  Array ft(shape); // top
  Array fb(shape); // bottom

  Array u(shape);
  Array v(shape);

  // continuous rainfall increment
  Array rain_step;
  if (p_rain_map && rain_rate > 0.f)
  {
    rain_step = (*p_rain_map) * (rain_rate * dt);
  }

  // --- instantiate runners once outside the iteration loop

  auto run_fp = clwrapper::Run("hydraulic_vpipes_flow_pass");

  run_fp.bind_imagef("z", z.vector, shape.x, shape.y); // inputs
  run_fp.bind_imagef("fl", fl.vector, shape.x, shape.y);
  run_fp.bind_imagef("fr", fr.vector, shape.x, shape.y);
  run_fp.bind_imagef("ft", ft.vector, shape.x, shape.y);
  run_fp.bind_imagef("fb", fb.vector, shape.x, shape.y);
  run_fp.bind_imagef("d1", d.vector, shape.x, shape.y);

  run_fp.bind_imagef("fl_out", fl.vector, shape.x, shape.y, true); // outputs
  run_fp.bind_imagef("fr_out", fr.vector, shape.x, shape.y, true);
  run_fp.bind_imagef("ft_out", ft.vector, shape.x, shape.y, true);
  run_fp.bind_imagef("fb_out", fb.vector, shape.x, shape.y, true);

  run_fp.bind_arguments(shape.x,
                        shape.y,
                        dt,
                        flux_diffusion ? 1 : 0,
                        flux_diffusion_strength,
                        outflow_boundaries ? 1 : 0);

  auto run_wa = clwrapper::Run("hydraulic_vpipes_water_pass");

  run_wa.bind_imagef("z", z.vector, shape.x, shape.y); // inputs
  run_wa.bind_imagef("fl", fl.vector, shape.x, shape.y);
  run_wa.bind_imagef("fr", fr.vector, shape.x, shape.y);
  run_wa.bind_imagef("ft", ft.vector, shape.x, shape.y);
  run_wa.bind_imagef("fb", fb.vector, shape.x, shape.y);
  run_wa.bind_imagef("d1", d.vector, shape.x, shape.y);

  run_wa.bind_imagef("d2_out", d.vector, shape.x, shape.y, true); // outputs
  run_wa.bind_imagef("u_out", u.vector, shape.x, shape.y, true);
  run_wa.bind_imagef("v_out", v.vector, shape.x, shape.y, true);

  run_wa.bind_arguments(shape.x,
                        shape.y,
                        dt,
                        water_height,
                        evap_rate,
                        outflow_boundaries ? 1 : 0);

  for (int it = 0; it < iterations; ++it)
  {
    // add continuous rainfall
    if (p_rain_map && rain_rate > 0.f)
    {
      d += rain_step;
    }
    else if (rain_rate > 0.f)
    {
      d += rain_rate * dt;
    }

    // --- flux update
    run_fp.write_imagef("fl");
    run_fp.write_imagef("fr");
    run_fp.write_imagef("ft");
    run_fp.write_imagef("fb");
    run_fp.write_imagef("d1");

    run_fp.execute({shape.x, shape.y});

    // update flux (from GPU to CPU)
    run_fp.read_imagef("fl_out");
    run_fp.read_imagef("fr_out");
    run_fp.read_imagef("ft_out");
    run_fp.read_imagef("fb_out");

    // --- water transport
    run_wa.write_imagef("fl");
    run_wa.write_imagef("fr");
    run_wa.write_imagef("ft");
    run_wa.write_imagef("fb");
    run_wa.write_imagef("d1");

    run_wa.execute({shape.x, shape.y});

    run_wa.read_imagef("d2_out");
  }

  // retrieve velocity field if requested
  if (p_vel_u)
  {
    run_wa.read_imagef("u_out");
    *p_vel_u = u;
  }
  if (p_vel_v)
  {
    run_wa.read_imagef("v_out");
    *p_vel_v = v;
  }

  // remove thin layer of remaining water
  if (dry_out_ratio != 0.f)
  {
    float dmax = d.max();
    water_depth_dry_out(d, dry_out_ratio, nullptr, dmax);
  }

  return d; // water depth
}

Array flow_simulation_viscous(const Array &z,
                              float        water_height,
                              const Array &depth_map,
                              int          iterations,
                              float        dt,
                              float        dry_out_ratio,
                              float        viscosity,
                              float        power,
                              float        evap_rate,
                              bool         outflow_boundaries)
{
  if (!validate_non_empty(z)) return Array();
  if (!validate_same_shape(z, depth_map)) return Array();

  const glm::ivec2 shape = z.shape;

  Array d = water_height * depth_map;

  auto run = clwrapper::Run("shallow_viscous_flow");

  run.bind_imagef("z", z.vector, shape.x, shape.y); // inputs
  run.bind_imagef("d_in", d.vector, shape.x, shape.y);
  run.bind_imagef("d_out", d.vector, shape.x, shape.y, true); // outputs

  run.bind_arguments(shape.x,
                     shape.y,
                     dt,
                     viscosity,
                     power,
                     evap_rate,
                     outflow_boundaries ? 1 : 0);

  run.write_imagef("z");

  for (int it = 0; it < iterations; ++it)
  {
    if (it % 10 == 0)
    {
      dt = helper_vflow_compute_adaptive_dt(d, viscosity, power);
      run.set_argument(5, dt);
    }

    run.write_imagef("d_in");

    run.execute({shape.x, shape.y});

    // update flux (from GPU to CPU)
    run.read_imagef("d_out");
  }

  // remove thin layer of remaining water
  if (dry_out_ratio != 0.f)
  {
    float dmax = d.max();
    water_depth_dry_out(d, dry_out_ratio, nullptr, dmax);
  }

  return d; // depth
}

} // namespace hmap::gpu
