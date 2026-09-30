R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel hydraulic_vpipes_flow_pass(read_only image2d_t  z,
                                       read_only image2d_t  fl,
                                       read_only image2d_t  fr,
                                       read_only image2d_t  ft,
                                       read_only image2d_t  fb,
                                       read_only image2d_t  d1,
                                       write_only image2d_t fl_out,
                                       write_only image2d_t fr_out,
                                       write_only image2d_t ft_out,
                                       write_only image2d_t fb_out,
                                       const int            nx,
                                       const int            ny,
                                       const float          dt,
                                       const int            flux_diffusion,
                                       const float flux_diffusion_strength,
                                       const int   outflow_boundaries)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const int   i = g.x;
  const int   j = g.y;
  const float plength = 1.f;
  const float gv = 1.f; // gravity

  float d0 = TGET(d1, i, j);
  float fl_val = TGET(fl, i, j);
  float fr_val = TGET(fr, i, j);
  float ft_val = TGET(ft, i, j);
  float fb_val = TGET(fb, i, j);

  // fast early-exit: if cell is dry and has no existing flux, outgoing fluxes
  // are strictly zero
  if (d0 <= 0.f && fl_val <= 0.f && fr_val <= 0.f && ft_val <= 0.f &&
      fb_val <= 0.f)
  {
    TSET(fl_out, i, j, 0.f);
    TSET(fr_out, i, j, 0.f);
    TSET(ft_out, i, j, 0.f);
    TSET(fb_out, i, j, 0.f);
    return;
  }

  float h0 = TGET(z, i, j) + d0;

  float dhl = (outflow_boundaries && i == 0)
                  ? d0
                  : (h0 - TGET(z, i - 1, j) - TGET(d1, i - 1, j));
  float dhr = (outflow_boundaries && i == nx - 1)
                  ? d0
                  : (h0 - TGET(z, i + 1, j) - TGET(d1, i + 1, j));
  float dht = (outflow_boundaries && j == ny - 1)
                  ? d0
                  : (h0 - TGET(z, i, j + 1) - TGET(d1, i, j + 1));
  float dhb = (outflow_boundaries && j == 0)
                  ? d0
                  : (h0 - TGET(z, i, j - 1) - TGET(d1, i, j - 1));

  float fl_new = max(0.f, fl_val + dt * gv * dhl / plength);
  float fr_new = max(0.f, fr_val + dt * gv * dhr / plength);
  float ft_new = max(0.f, ft_val + dt * gv * dht / plength);
  float fb_new = max(0.f, fb_val + dt * gv * dhb / plength);

  // flux diffusion
  if (flux_diffusion)
  {
    float c0 = flux_diffusion_strength;
    float c1 = 1.f - 4.f * c0;
    float f_diff = c0 * (fl_new + fr_new + ft_new + fb_new);

    fl_new = c1 * fl_new + f_diff;
    fr_new = c1 * fr_new + f_diff;
    ft_new = c1 * ft_new + f_diff;
    fb_new = c1 * fb_new + f_diff;
  }

  // normalize
  float sum = fl_new + fr_new + ft_new + fb_new;
  float k = 0.f;
  if (sum > 1e-5f) k = d0 * plength * plength / (sum * dt);
  k = clamp(k, 0.f, 1.f);

  // output
  TSET(fl_out, i, j, fl_new * k);
  TSET(fr_out, i, j, fr_new * k);
  TSET(ft_out, i, j, ft_new * k);
  TSET(fb_out, i, j, fb_new * k);
}

void kernel hydraulic_vpipes_water_pass(read_only image2d_t  z,
                                        read_only image2d_t  fl,
                                        read_only image2d_t  fr,
                                        read_only image2d_t  ft,
                                        read_only image2d_t  fb,
                                        read_only image2d_t  d1,
                                        write_only image2d_t d2_out,
                                        write_only image2d_t u_out,
                                        write_only image2d_t v_out,
                                        const int            nx,
                                        const int            ny,
                                        const float          dt,
                                        const float          water_height,
                                        const float          evap_rate,
                                        const int            outflow_boundaries)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const int   i = g.x;
  const int   j = g.y;
  const float plength = 1.f;
  const float gv = 1.f; // gravity

  float d0 = TGET(d1, i, j);
  float fl_val = TGET(fl, i, j);
  float fr_val = TGET(fr, i, j);
  float ft_val = TGET(ft, i, j);
  float fb_val = TGET(fb, i, j);

  float in_l = (outflow_boundaries && i == 0) ? 0.f : TGET(fr, i - 1, j);
  float in_r = (outflow_boundaries && i == nx - 1) ? 0.f : TGET(fl, i + 1, j);
  float in_b = (outflow_boundaries && j == 0) ? 0.f : TGET(ft, i, j - 1);
  float in_t = (outflow_boundaries && j == ny - 1) ? 0.f : TGET(fb, i, j + 1);

  // fast early-exit: if cell is dry and has no incoming or outgoing fluxes
  if (d0 <= 0.f && fl_val <= 0.f && fr_val <= 0.f && ft_val <= 0.f &&
      fb_val <= 0.f && in_l <= 0.f && in_r <= 0.f && in_b <= 0.f && in_t <= 0.f)
  {
    TSET(d2_out, i, j, 0.f);
    TSET(u_out, i, j, 0.f);
    TSET(v_out, i, j, 0.f);
    return;
  }

  float dv = dt *
             (in_l + in_b + in_r + in_t - fl_val - fr_val - ft_val - fb_val);
  float d2_new = max(0.f, d0 + dv / (plength * plength));

  if (evap_rate > 0.f)
  {
    d2_new = max(0.f, d2_new * (1.f - evap_rate * dt));
  }

  TSET(d2_out, i, j, d2_new);

  float u_new = 0.5f *
                (in_l - fl_val + fr_val -
                 ((outflow_boundaries && i == nx - 1) ? 0.f
                                                      : TGET(fl, i + 1, j)));
  float v_new = 0.5f *
                (in_b - fb_val + ft_val -
                 ((outflow_boundaries && j == ny - 1) ? 0.f
                                                      : TGET(fb, i, j + 1)));

  float dmean = max(0.001f * water_height, d2_new);

  u_new /= dmean;
  v_new /= dmean;

  TSET(u_out, i, j, u_new);
  TSET(v_out, i, j, v_new);
}

void kernel
hydraulic_vpipes_erosion_pass(read_only image2d_t  z,
                              read_only image2d_t  d2,
                              read_only image2d_t  u,
                              read_only image2d_t  v,
                              read_only image2d_t  s,
                              write_only image2d_t z_out,
                              write_only image2d_t s_out,
                              const int            nx,
                              const int            ny,
                              const float          water_height,
                              const float          k_capacity,
                              const float          k_erode,
                              const float          k_depose,
                              const float          k_discharge_exp,
                              const float          downcutting_max_depth_ratio)
{
  const float salpha_min = 0.01f;

  //

  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const int i = g.x;
  const int j = g.y;

  // gradient norm
  float dzx = 0.5f * (TGET(z, i + 1, j) - TGET(z, i - 1, j) +
                      TGET(s, i + 1, j) - TGET(s, i - 1, j));
  float dzy = 0.5f * (TGET(z, i, j + 1) - TGET(z, i, j - 1) +
                      TGET(s, i, j + 1) - TGET(s, i, j - 1));

  float talus = hypot(dzx, dzy);
  float dzn = max(1e-3f, nx * talus);

  // sediment capacity
  const float dmax = water_height;

  float salpha = max(salpha_min, dzn / hypot(1.f, dzn));
  float speed = hypot(TGET(u, i, j), TGET(v, i, j));
  float depth = min(TGET(d2, i, j) / dmax, downcutting_max_depth_ratio);
  float discharge = depth * speed;
  float capa = k_capacity * pow(discharge, k_discharge_exp) * salpha;

  float st = TGET(s, i, j);
  float zt = TGET(z, i, j);

  if (capa > st)
  {
    // erosion
    float amount = k_erode * (capa - st);
    TSET(z_out, i, j, zt - amount);
    TSET(s_out, i, j, st + amount);
  }
  else
  {
    // deposition
    float amount = k_depose * (st - capa);
    TSET(z_out, i, j, zt + amount);
    TSET(s_out, i, j, st - amount);
  }
}

void kernel hydraulic_vpipes_sediment_transport_pass(read_only image2d_t  u,
                                                     read_only image2d_t  v,
                                                     read_only image2d_t  s,
                                                     write_only image2d_t s_out,
                                                     const int            nx,
                                                     const int            ny,
                                                     const float          dt)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const sampler_t sampler_itp = CLK_NORMALIZED_COORDS_FALSE |
                                CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  const int i = g.x;
  const int j = g.y;

  float2 pos = (float2)(i + 0.5f - dt * TGET(u, i, j),
                        j + 0.5f - dt * TGET(v, i, j));
  float  s_new = read_imagef(s, sampler_itp, pos).x;

  TSET(s_out, i, j, s_new);
}

// Adds one rainfall increment. FP contraction is disabled so the result is
// bit-identical to the former host-side `d += rain_map * amount`.
#pragma OPENCL FP_CONTRACT OFF
void kernel hydraulic_vpipes_rain_pass(read_only image2d_t  d_in,
                                       read_only image2d_t  rain_map,
                                       write_only image2d_t d_out,
                                       const int            nx,
                                       const int            ny,
                                       const float          amount,
                                       const int            use_map)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float d = TGET(d_in, g.x, g.y);
  float r = use_map ? TGET(rain_map, g.x, g.y) * amount : amount;

  TSET(d_out, g.x, g.y, d + r);
}
#pragma OPENCL FP_CONTRACT DEFAULT
)""
