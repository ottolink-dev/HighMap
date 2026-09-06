R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

kernel void shallow_viscous_flow(read_only image2d_t  z,
                                 read_only image2d_t  h_in,
                                 write_only image2d_t h_out,
                                 int                  nx,
                                 int                  ny,
                                 float                dt,
                                 float                viscosity,
                                 float                power,
                                 float                evap_rate,
                                 int                  outflow_boundaries)
{
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i >= nx || j >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
  const float diag = 0.70710678f;

  float h = TGET(h_in, i, j);
  float hxp = (outflow_boundaries && i == nx - 1) ? 0.f : TGET(h_in, i + 1, j);
  float hxm = (outflow_boundaries && i == 0) ? 0.f : TGET(h_in, i - 1, j);
  float hyp = (outflow_boundaries && j == ny - 1) ? 0.f : TGET(h_in, i, j + 1);
  float hym = (outflow_boundaries && j == 0) ? 0.f : TGET(h_in, i, j - 1);

  // fast early-exit: skip all elevation fetches, pow(), and flux math if cell
  // and neighbors are dry
  if (h <= 0.f && hxp <= 0.f && hxm <= 0.f && hyp <= 0.f && hym <= 0.f)
  {
    TSET(h_out, i, j, 0.f);
    return;
  }

  float zc = TGET(z, i, j);

  float Hc = h + zc;
  float Hxp = (outflow_boundaries && i == nx - 1) ? zc
                                                  : (hxp + TGET(z, i + 1, j));
  float Hxm = (outflow_boundaries && i == 0) ? zc : (hxm + TGET(z, i - 1, j));
  float Hyp = (outflow_boundaries && j == ny - 1) ? zc
                                                  : (hyp + TGET(z, i, j + 1));
  float Hym = (outflow_boundaries && j == 0) ? zc : (hym + TGET(z, i, j - 1));

  // upwind face mobilities
  float h_up_xp = (Hc >= Hxp) ? h : hxp;
  float h_up_xm = (Hxm >= Hc) ? hxm : h;
  float h_up_yp = (Hc >= Hyp) ? h : hyp;
  float h_up_ym = (Hym >= Hc) ? hym : h;

  float Mp = pow(h_up_xp, power) / viscosity;
  float Mm = pow(h_up_xm, power) / viscosity;
  float Myp = pow(h_up_yp, power) / viscosity;
  float Mym = pow(h_up_ym, power) / viscosity;

  // fluxes through faces
  float qxp = -Mp * (Hxp - Hc);
  float qxm = -Mm * (Hc - Hxm);
  float qyp = -Myp * (Hyp - Hc);
  float qym = -Mym * (Hc - Hym);

  // outflow limiting to preserve non-negativity and mass conservation
  float q_out = max(0.f, qxp) + max(0.f, -qxm) + max(0.f, qyp) + max(0.f, -qym);
  float q_in = max(0.f, -qxp) + max(0.f, qxm) + max(0.f, -qyp) + max(0.f, qym);

  float k_out = 1.f;
  if (q_out * dt > h && q_out > 1e-6f)
  {
    k_out = h / (q_out * dt);
  }

  // explicit update
  float h_new = max(0.f, h + dt * (q_in - k_out * q_out));

  if (evap_rate > 0.f)
  {
    h_new = max(0.f, h_new * (1.f - evap_rate * dt));
  }

  // smoothing
  float h_avg = (TGET(h_in, i - 1, j) + TGET(h_in, i + 1, j) +
                 TGET(h_in, i, j - 1) + TGET(h_in, i, j + 1) +
                 diag * (TGET(h_in, i - 1, j - 1) + TGET(h_in, i + 1, j - 1) +
                         TGET(h_in, i - 1, j + 1) + TGET(h_in, i + 1, j + 1))) /
                (4.f + 4.f * diag);

  float k_visc = 0.001f;
  // h_new = mix(h_new, h_avg, k_visc);

  TSET(h_out, i, j, h_new);
}
)""
