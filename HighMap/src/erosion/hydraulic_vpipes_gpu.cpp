/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/internal/validation.hpp"

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

  Array rain_map(shape, 1.f);
  if (p_rain_map) rain_map = *p_rain_map;

  Array d(shape, water_height); // water height
  Array d1(shape);
  Array d2(shape);

  Array s(shape); // sediment height

  Array fl(shape); // left flux
  Array fr(shape); // right
  Array ft(shape); // top
  Array fb(shape); // bottom

  Array u(shape);
  Array v(shape);

  d *= rain_map;
  const float water_volume_init = d.sum();
  const float rain_map_volume = rain_map.sum();

  // --- Main loop

  for (int it = 0; it < iterations; ++it)
  {
    // water volume increment
    if (maintain_water_volume)
    {
      float water_volume = d.sum();
      float rain_rate = (water_volume_init - water_volume) / rain_map_volume;
      d += rain_rate * rain_map;
    }

    d1 = d;

    // --- Flux update

    auto run_fp = clwrapper::Run("hydraulic_vpipes_flow_pass");

    run_fp.bind_imagef("z", z.vector, shape.x, shape.y);
    run_fp.bind_imagef("fl", fl.vector, shape.x, shape.y);
    run_fp.bind_imagef("fr", fr.vector, shape.x, shape.y);
    run_fp.bind_imagef("ft", ft.vector, shape.x, shape.y);
    run_fp.bind_imagef("fb", fb.vector, shape.x, shape.y);
    run_fp.bind_imagef("d1", d1.vector, shape.x, shape.y);

    run_fp.bind_imagef("fl_out", fl.vector, shape.x, shape.y, true);
    run_fp.bind_imagef("fr_out", fr.vector, shape.x, shape.y, true);
    run_fp.bind_imagef("ft_out", ft.vector, shape.x, shape.y, true);
    run_fp.bind_imagef("fb_out", fb.vector, shape.x, shape.y, true);

    run_fp.bind_arguments(shape.x,
                          shape.y,
                          dt,
                          flux_diffusion ? 1 : 0,
                          flux_diffusion_strength,
                          0);

    run_fp.execute({shape.x, shape.y});

    run_fp.read_imagef("fl_out");
    run_fp.read_imagef("fr_out");
    run_fp.read_imagef("ft_out");
    run_fp.read_imagef("fb_out");

    // --- Water transport

    auto run_wa = clwrapper::Run("hydraulic_vpipes_water_pass");

    run_wa.bind_imagef("z", z.vector, shape.x, shape.y);
    run_wa.bind_imagef("fl", fl.vector, shape.x, shape.y);
    run_wa.bind_imagef("fr", fr.vector, shape.x, shape.y);
    run_wa.bind_imagef("ft", ft.vector, shape.x, shape.y);
    run_wa.bind_imagef("fb", fb.vector, shape.x, shape.y);
    run_wa.bind_imagef("d1", d1.vector, shape.x, shape.y);

    run_wa.bind_imagef("d2_out", d2.vector, shape.x, shape.y, true);
    run_wa.bind_imagef("u_out", u.vector, shape.x, shape.y, true);
    run_wa.bind_imagef("v_out", v.vector, shape.x, shape.y, true);

    run_wa.bind_arguments(shape.x, shape.y, dt, water_height, evap_rate, 0);

    run_wa.execute({shape.x, shape.y});

    run_wa.read_imagef("d2_out");
    run_wa.read_imagef("u_out");
    run_wa.read_imagef("v_out");

    // --- Erosion and deposition

    auto run_er = clwrapper::Run("hydraulic_vpipes_erosion_pass");

    run_er.bind_imagef("z", z.vector, shape.x, shape.y);
    run_er.bind_imagef("d2", d2.vector, shape.x, shape.y);
    run_er.bind_imagef("u", u.vector, shape.x, shape.y);
    run_er.bind_imagef("v", v.vector, shape.x, shape.y);
    run_er.bind_imagef("s", s.vector, shape.x, shape.y);

    run_er.bind_imagef("z_out", z.vector, shape.x, shape.y, true);
    run_er.bind_imagef("s_out", s.vector, shape.x, shape.y, true);

    run_er.bind_arguments(shape.x,
                          shape.y,
                          water_height,
                          k_capacity,
                          k_erode,
                          k_depose,
                          k_discharge_exp,
                          downcutting_max_depth_ratio);

    run_er.execute({shape.x, shape.y});

    run_er.read_imagef("z_out");
    run_er.read_imagef("s_out");

    // --- Sediment transport

    auto run_st = clwrapper::Run("hydraulic_vpipes_sediment_transport_pass");

    run_st.bind_imagef("u", u.vector, shape.x, shape.y);
    run_st.bind_imagef("v", v.vector, shape.x, shape.y);
    run_st.bind_imagef("s", s.vector, shape.x, shape.y);

    run_st.bind_imagef("s_out", s.vector, shape.x, shape.y, true);

    run_st.bind_arguments(shape.x, shape.y, dt);

    run_st.execute({shape.x, shape.y});

    run_st.read_imagef("s_out");

    // update state variable
    d = d2;
  }

  // --- Outputs

  if (p_water_depth) *p_water_depth = d;
  if (p_sediment) *p_sediment = s;
  if (p_vel_u) *p_vel_u = u;
  if (p_vel_v) *p_vel_v = v;
}

} // namespace hmap::gpu
