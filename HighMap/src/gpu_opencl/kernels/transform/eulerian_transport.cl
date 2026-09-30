R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel eulerian_transport(read_only image2d_t  u,
                               read_only image2d_t  v,
                               read_only image2d_t  field_in,
                               write_only image2d_t field_out,
                               const int            nx,
                               const int            ny,
                               const float          dt)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  int  idx = linear_index(g.x, g.y, nx);

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
  const int i = g.x;
  const int j = g.y;

  float A = TGET(field_in, i, j);
  float u0 = TGET(u, i, j);
  float v0 = TGET(v, i, j);

  // neighbor values
  float AL = TGET(field_in, i - 1, j);
  float AR = TGET(field_in, i + 1, j);
  float AD = TGET(field_in, i, j - 1);
  float AU = TGET(field_in, i, j + 1);

  float uL = TGET(u, i - 1, j);
  float uR = TGET(u, i + 1, j);
  float vD = TGET(v, i, j - 1);
  float vU = TGET(v, i, j + 1);

  // fluxes
  float FxL = AL * uL;
  float FxR = AR * uR;
  float FyD = AD * vD;
  float FyU = AU * vU;

  // divergence of flux
  float div = (FxR - FxL + FyU - FyD) * 0.5f; // central difference

  // upwind scheme for stabilization
  float inflow = 0.f;
  float outflow = 0.f;

  if (u0 > 0)
  {
    inflow += TGET(field_in, i - 1, j) * u0;
    outflow += A * u0;
  }
  else
  {
    inflow += TGET(field_in, i + 1, j) * (-u0);
    outflow += A * (-u0);
  }

  if (v0 > 0)
  {
    inflow += TGET(field_in, i, j - 1) * v0;
    outflow += A * v0;
  }
  else
  {
    inflow += TGET(field_in, i, j + 1) * (-v0);
    outflow += A * (-v0);
  }

  float src = hypot(u0, v0);
  float A_new = A + dt * (src + inflow - outflow);

  /* // source term (velocity-only field_inumulation) */
  /* float src = hypot(u0, v0); */

  /* // update */
  /* float A_new = A - dt * div + dt * src; */

  TSET(field_out, i, j, A_new);
}
)""
