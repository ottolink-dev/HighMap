/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/erosion/hydraulic_erosion.hpp"
#include "highmap/internal/validation.hpp"

#include <CL/opencl.hpp>

namespace hmap::gpu
{

void hydraulic_vpipes(Array &z,
                      float  water_height,
                      bool   maintain_water_volume,
                      float  evap_rate,
                      int    iterations,
                      float  dt,
                      float  k_capacity,
                      float  k_erode,
                      float  k_depose,
                      float  k_discharge_exp,
                      float  downcutting_max_depth_ratio,
                      bool   flux_diffusion,
                      float  flux_diffusion_strength,
                      Array *p_rain_map,
                      Array *p_water_depth,
                      Array *p_sediment,
                      Array *p_vel_u,
                      Array *p_vel_v)
{
  if (!validate_non_empty(z)) return;
  if (p_rain_map && !validate_same_shape(z, *p_rain_map)) return;

  const glm::ivec2 shape = z.shape;

  if (iterations <= 0)
  {
    if (p_water_depth) *p_water_depth = Array(shape, 0.f);
    if (p_sediment) *p_sediment = Array(shape, 0.f);
    if (p_vel_u) *p_vel_u = Array(shape, 0.f);
    if (p_vel_v) *p_vel_v = Array(shape, 0.f);
    return;
  }

  // --- Device state: elevation, depth, sediment, and flux ping-pong pairs

  Array zeros(shape, 0.f);
  Array d = water_height * (p_rain_map ? *p_rain_map : Array(shape, 1.f));
  Array s(shape, 0.f);
  Array u(shape);
  Array v(shape);

  using clwrapper::Direction;

  // flux pass: (z, fl, fr, ft, fb, d1, fl_out, fr_out, ft_out, fb_out, ...)
  auto run_fp = clwrapper::Run("hydraulic_vpipes_flow_pass");

  run_fp.bind_imagef("z_a", z.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("fl_a", zeros.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("fr_a", zeros.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("ft_a", zeros.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("fb_a", zeros.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("d_a", d.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("fl_b", zeros.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("fr_b", zeros.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("ft_b", zeros.vector, shape.x, shape.y, Direction::INOUT);
  run_fp.bind_imagef("fb_b", zeros.vector, shape.x, shape.y, Direction::INOUT);

  run_fp.bind_arguments(shape.x,
                        shape.y,
                        dt,
                        flux_diffusion ? 1 : 0,
                        flux_diffusion_strength,
                        0);

  // water pass: (z, fl, fr, ft, fb, d1, d2_out, u_out, v_out, ...)
  auto run_wa = clwrapper::Run("hydraulic_vpipes_water_pass",
                               run_fp.get_queue());

  run_wa.bind_image2d("z_a", run_fp.get_image2d("z_a"));
  run_wa.bind_image2d("fl_b", run_fp.get_image2d("fl_b"));
  run_wa.bind_image2d("fr_b", run_fp.get_image2d("fr_b"));
  run_wa.bind_image2d("ft_b", run_fp.get_image2d("ft_b"));
  run_wa.bind_image2d("fb_b", run_fp.get_image2d("fb_b"));
  run_wa.bind_image2d("d_a", run_fp.get_image2d("d_a"));
  run_wa.bind_imagef("d_b", d.vector, shape.x, shape.y, Direction::INOUT);
  run_wa.bind_imagef("u", u.vector, shape.x, shape.y, Direction::INOUT);
  run_wa.bind_imagef("v", v.vector, shape.x, shape.y, Direction::INOUT);

  run_wa.bind_arguments(shape.x, shape.y, dt, water_height, evap_rate, 0);

  // erosion pass: (z, d2, u, v, s, z_out, s_out, ...)
  auto run_er = clwrapper::Run("hydraulic_vpipes_erosion_pass",
                               run_fp.get_queue());

  run_er.bind_image2d("z_a", run_fp.get_image2d("z_a"));
  run_er.bind_image2d("d_b", run_wa.get_image2d("d_b"));
  run_er.bind_image2d("u", run_wa.get_image2d("u"));
  run_er.bind_image2d("v", run_wa.get_image2d("v"));
  run_er.bind_imagef("s_a", s.vector, shape.x, shape.y, Direction::INOUT);
  run_er.bind_imagef("z_b", z.vector, shape.x, shape.y, Direction::INOUT);
  run_er.bind_imagef("s_b", s.vector, shape.x, shape.y, Direction::INOUT);

  run_er.bind_arguments(shape.x,
                        shape.y,
                        water_height,
                        k_capacity,
                        k_erode,
                        k_depose,
                        k_discharge_exp,
                        downcutting_max_depth_ratio);

  // sediment transport pass: (u, v, s, s_out, ...)
  auto run_st = clwrapper::Run("hydraulic_vpipes_sediment_transport_pass",
                               run_fp.get_queue());

  run_st.bind_image2d("u", run_wa.get_image2d("u"));
  run_st.bind_image2d("v", run_wa.get_image2d("v"));
  run_st.bind_image2d("s_b", run_er.get_image2d("s_b"));
  run_st.bind_image2d("s_a", run_er.get_image2d("s_a"));

  run_st.bind_arguments(shape.x, shape.y, dt);

  // rain pass (optional, if maintain_water_volume && evap_rate > 0.f)
  const bool use_rain = maintain_water_volume && evap_rate > 0.f;
  const bool use_map = use_rain && (p_rain_map != nullptr);

  std::unique_ptr<clwrapper::Run> run_rain;
  if (use_rain)
  {
    run_rain = std::make_unique<clwrapper::Run>("hydraulic_vpipes_rain_pass",
                                                run_fp.get_queue());

    run_rain->bind_image2d("d_a", run_fp.get_image2d("d_a"));
    run_rain->bind_imagef("rain",
                          use_map ? p_rain_map->vector : zeros.vector,
                          shape.x,
                          shape.y);
    run_rain->bind_image2d("d_b", run_wa.get_image2d("d_b"));
    run_rain->bind_arguments(shape.x,
                             shape.y,
                             water_height * evap_rate * dt,
                             use_map ? 1 : 0);
  }

  // ping-pong handles
  const std::array<cl::Image2D, 2> img_z = {run_fp.get_image2d("z_a").cl_image,
                                            run_er.get_image2d("z_b").cl_image};

  const std::array<cl::Image2D, 2> img_d = {run_fp.get_image2d("d_a").cl_image,
                                            run_wa.get_image2d("d_b").cl_image};

  const std::array<cl::Image2D, 2> img_s = {run_er.get_image2d("s_a").cl_image,
                                            run_er.get_image2d("s_b").cl_image};

  const std::array<std::array<cl::Image2D, 4>, 2> img_f = {
      {{run_fp.get_image2d("fl_a").cl_image,
        run_fp.get_image2d("fr_a").cl_image,
        run_fp.get_image2d("ft_a").cl_image,
        run_fp.get_image2d("fb_a").cl_image},
       {run_fp.get_image2d("fl_b").cl_image,
        run_fp.get_image2d("fr_b").cl_image,
        run_fp.get_image2d("ft_b").cl_image,
        run_fp.get_image2d("fb_b").cl_image}}};

  int zc = 0; // index of current elevation image
  int dc = 0; // index of current depth image
  int fc = 0; // index of current flux images
  int sc = 0; // index of current sediment image

  // --- Main loop

  for (int it = 0; it < iterations; ++it)
  {
    // continuous rainfall: d[dc] + rain -> d[1 - dc]
    if (use_rain)
    {
      run_rain->set_argument(0, img_d[dc]);
      run_rain->set_argument(2, img_d[1 - dc]);
      run_rain->execute_async({shape.x, shape.y});
      dc = 1 - dc;
    }

    // flux update: reads z[zc], f[fc], d[dc]; writes f[1 - fc]
    run_fp.set_argument(0, img_z[zc]);
    for (int k = 0; k < 4; ++k)
      run_fp.set_argument(1 + k, img_f[fc][k]);
    run_fp.set_argument(5, img_d[dc]);
    for (int k = 0; k < 4; ++k)
      run_fp.set_argument(6 + k, img_f[1 - fc][k]);

    run_fp.execute_async({shape.x, shape.y});

    // water transport: reads z[zc], f[1 - fc], d[dc]; writes d[1 - dc], u, v
    run_wa.set_argument(0, img_z[zc]);
    for (int k = 0; k < 4; ++k)
      run_wa.set_argument(1 + k, img_f[1 - fc][k]);
    run_wa.set_argument(5, img_d[dc]);
    run_wa.set_argument(6, img_d[1 - dc]);

    run_wa.execute_async({shape.x, shape.y});

    fc = 1 - fc;
    dc = 1 - dc;

    // erosion and deposition: reads z[zc], d[dc], u, v, s[sc]; writes z[1 -
    // zc], s[1 - sc]
    run_er.set_argument(0, img_z[zc]);
    run_er.set_argument(1, img_d[dc]);
    run_er.set_argument(4, img_s[sc]);
    run_er.set_argument(5, img_z[1 - zc]);
    run_er.set_argument(6, img_s[1 - sc]);

    run_er.execute_async({shape.x, shape.y});

    zc = 1 - zc;
    sc = 1 - sc;

    // sediment transport: reads u, v, s[sc]; writes s[1 - sc]
    run_st.set_argument(2, img_s[sc]);
    run_st.set_argument(3, img_s[1 - sc]);

    run_st.execute_async({shape.x, shape.y});

    sc = 1 - sc;
  }

  run_st.finish();

  // --- Outputs

  if (zc == 0)
    run_fp.read_imagef("z_a");
  else
    run_er.read_imagef("z_b");

  if (p_water_depth)
  {
    if (dc == 0)
      run_fp.read_imagef("d_a");
    else
      run_wa.read_imagef("d_b");
    *p_water_depth = d;
  }

  if (p_sediment)
  {
    if (sc == 0)
      run_er.read_imagef("s_a");
    else
      run_er.read_imagef("s_b");
    *p_sediment = s;
  }

  if (p_vel_u)
  {
    run_wa.read_imagef("u");
    *p_vel_u = u;
  }

  if (p_vel_v)
  {
    run_wa.read_imagef("v");
    *p_vel_v = v;
  }
}

} // namespace hmap::gpu
