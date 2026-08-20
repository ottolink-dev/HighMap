/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <cmath>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/hydrology/hydrology.hpp"

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
                      float        dry_out_ratio)
{
  const glm::ivec2 shape = z.shape;

  Array d = water_height * depth_map;
  Array fl(shape); // left flux
  Array fr(shape); // right
  Array ft(shape); // top
  Array fb(shape); // bottom

  Array u(shape); // velocity (output of the water pass)
  Array v(shape);

  // All the state (terrain, water depth, fluxes) is uploaded once and then
  // kept GPU-resident for the whole simulation using ping-pong images; the
  // host is only touched again for the final readback.

  clwrapper::Run run_fp("hydraulic_vpipes_flow_pass");

  // inputs (slot 0-5)
  run_fp.bind_imagef("z", z.vector, shape.x, shape.y);
  run_fp.bind_imagef("fl",
                     fl.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);
  run_fp.bind_imagef("fr",
                     fr.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);
  run_fp.bind_imagef("ft",
                     ft.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);
  run_fp.bind_imagef("fb",
                     fb.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);
  run_fp.bind_imagef("d1",
                     d.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);

  // outputs (slot 6-9), ping-pong counterparts of the flux inputs
  run_fp.bind_imagef("fl_out",
                     fl.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);
  run_fp.bind_imagef("fr_out",
                     fr.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);
  run_fp.bind_imagef("ft_out",
                     ft.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);
  run_fp.bind_imagef("fb_out",
                     fb.vector,
                     shape.x,
                     shape.y,
                     clwrapper::Direction::INOUT);

  // ping-pong counterpart of the water depth input (not a kernel argument
  // of this pass, only shared with the water pass below)
  cl::Image2D img_d_alt = run_fp.create_imagef("d2",
                                               d.vector,
                                               shape.x,
                                               shape.y,
                                               clwrapper::Direction::INOUT);

  run_fp.bind_arguments(shape.x,
                        shape.y,
                        dt,
                        flux_diffusion ? 1 : 0,
                        flux_diffusion_strength);

  // --- water transport, reads the device images of the flow pass directly

  clwrapper::Run run_wa("hydraulic_vpipes_water_pass");

  run_wa.bind_imagef("z", z.vector, shape.x, shape.y, run_fp.get_imagef("z"));
  run_wa.bind_imagef("fl",
                     fl.vector,
                     shape.x,
                     shape.y,
                     run_fp.get_imagef("fl_out"));
  run_wa.bind_imagef("fr",
                     fr.vector,
                     shape.x,
                     shape.y,
                     run_fp.get_imagef("fr_out"));
  run_wa.bind_imagef("ft",
                     ft.vector,
                     shape.x,
                     shape.y,
                     run_fp.get_imagef("ft_out"));
  run_wa.bind_imagef("fb",
                     fb.vector,
                     shape.x,
                     shape.y,
                     run_fp.get_imagef("fb_out"));
  run_wa.bind_imagef("d1",
                     d.vector,
                     shape.x,
                     shape.y,
                     run_fp.get_imagef("d1"));
  run_wa.bind_imagef("d2_out", d.vector, shape.x, shape.y, img_d_alt);
  run_wa.bind_imagef("u_out", u.vector, shape.x, shape.y, true); // outputs
  run_wa.bind_imagef("v_out", v.vector, shape.x, shape.y, true);

  run_wa.bind_arguments(shape.x, shape.y, dt, water_height);

  // --- main loop, swap image arguments to ping-pong between the two states

  const cl::Image2D img_flux_a[4] = {run_fp.get_imagef("fl"),
                                     run_fp.get_imagef("fr"),
                                     run_fp.get_imagef("ft"),
                                     run_fp.get_imagef("fb")};
  const cl::Image2D img_flux_b[4] = {run_fp.get_imagef("fl_out"),
                                     run_fp.get_imagef("fr_out"),
                                     run_fp.get_imagef("ft_out"),
                                     run_fp.get_imagef("fb_out")};
  const cl::Image2D img_d_main = run_fp.get_imagef("d1");

  for (int it = 0; it < iterations; ++it)
  {
    const int in_side = it & 1;
    const int out_side = in_side ^ 1;

    // flux update: reads state #in_side, writes state #out_side
    for (int k = 0; k < 4; ++k)
      run_fp.set_argument(1 + k,
                          in_side == 0 ? img_flux_a[k] : img_flux_b[k]);
    run_fp.set_argument(5, in_side == 0 ? img_d_main : img_d_alt);
    for (int k = 0; k < 4; ++k)
      run_fp.set_argument(6 + k,
                          out_side == 0 ? img_flux_a[k] : img_flux_b[k]);

    run_fp.execute({shape.x, shape.y});

    // water transport: consumes the fluxes just computed
    for (int k = 0; k < 4; ++k)
      run_wa.set_argument(1 + k,
                          out_side == 0 ? img_flux_a[k] : img_flux_b[k]);
    run_wa.set_argument(5, in_side == 0 ? img_d_main : img_d_alt);
    run_wa.set_argument(6, out_side == 0 ? img_d_main : img_d_alt);

    run_wa.execute({shape.x, shape.y});
  }

  // final readback of the water depth (the fluxes and velocities are never
  // needed on the host)
  if (iterations > 0)
  {
    if (iterations % 2 == 0)
      run_fp.read_imagef("d1");
    else
      run_fp.read_imagef("d2");
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
                              float        power)
{
  const glm::ivec2 shape = z.shape;

  Array d = water_height * depth_map;
  Array d_alt(shape, 0.f); // readback buffer of the ping-pong counterpart

  // The simulation state stays on the GPU: two depth images are used in
  // alternance (ping-pong), so one iteration is a single kernel launch with
  // zero host/device transfer. The terrain is uploaded only once.

  clwrapper::Run run_a("shallow_viscous_flow");

  run_a.bind_imagef("z", z.vector, shape.x, shape.y); // input
  run_a.bind_imagef("d_in",
                    d.vector,
                    shape.x,
                    shape.y,
                    clwrapper::Direction::INOUT);
  run_a.bind_imagef("d_out",
                    d_alt.vector,
                    shape.x,
                    shape.y,
                    clwrapper::Direction::INOUT);

  run_a.bind_arguments(shape.x, shape.y, dt, viscosity, power);

  clwrapper::Run run_b("shallow_viscous_flow");

  run_b.bind_imagef("z", z.vector, shape.x, shape.y); // input
  run_b.bind_imagef("d_in",
                    d_alt.vector,
                    shape.x,
                    shape.y,
                    run_a.get_imagef("d_out")); // cross-wired ping-pong
  run_b.bind_imagef("d_out",
                    d.vector,
                    shape.x,
                    shape.y,
                    run_a.get_imagef("d_in"));

  run_b.bind_arguments(shape.x, shape.y, dt, viscosity, power);

  Array       *p_d_last = &d;       // host mirror of the latest device state
  clwrapper::Run *p_run_last = nullptr; // who wrote it

  for (int it = 0; it < iterations; ++it)
  {
    // refresh the time step regularly (needs one readback of the current
    // depth map)
    if (it % 10 == 0)
    {
      if (p_run_last) p_run_last->read_imagef("d_out");

      dt = helper_vflow_compute_adaptive_dt(*p_d_last, viscosity, power);
      run_a.set_argument(5, dt);
      run_b.set_argument(5, dt);
    }

    if (it % 2 == 0)
    {
      run_a.execute({shape.x, shape.y});
      p_run_last = &run_a;
      p_d_last = &d_alt;
    }
    else
    {
      run_b.execute({shape.x, shape.y});
      p_run_last = &run_b;
      p_d_last = &d;
    }
  }

  // final readback
  if (p_run_last)
  {
    p_run_last->read_imagef("d_out");
    if (p_d_last != &d) d = std::move(d_alt);
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
