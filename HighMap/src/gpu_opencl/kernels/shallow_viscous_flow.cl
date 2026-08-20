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
                                 float                power)
{
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i >= nx || j >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float h = TGET(h_in, i, j);

  // neighbors
  float hxp = TGET(h_in, i + 1, j);
  float hxm = TGET(h_in, i - 1, j);
  float hyp = TGET(h_in, i, j + 1);
  float hym = TGET(h_in, i, j - 1);

  // a dry cell surrounded by dry cells cannot evolve (common case once the
  // fluid front has stabilized), skip the expensive mobility evaluation
  if (h == 0.f && hxp == 0.f && hxm == 0.f && hyp == 0.f && hym == 0.f)
  {
    TSET(h_out, i, j, 0.f);
    return;
  }

  float zc = TGET(z, i, j);

  // total water surface elevations
  float Hc = h + zc;
  float Hxp = hxp + TGET(z, i + 1, j);
  float Hxm = hxm + TGET(z, i - 1, j);
  float Hyp = hyp + TGET(z, i, j + 1);
  float Hym = hym + TGET(z, i, j - 1);

  // face mobilities (power-law fluid, viscosity = 1 -> water)
  float inv_visc = 1.f / viscosity;
  float Mp = pow(fmax(0.f, 0.5f * (h + hxp)), power) * inv_visc;
  float Mm = pow(fmax(0.f, 0.5f * (h + hxm)), power) * inv_visc;
  float Myp = pow(fmax(0.f, 0.5f * (h + hyp)), power) * inv_visc;
  float Mym = pow(fmax(0.f, 0.5f * (h + hym)), power) * inv_visc;

  // Outgoing flux through each face, capped so that the source cell cannot
  // give away more than a quarter of its depth per face and per time step.
  // Both cells sharing a face compute the same capped value (same mobility,
  // same elevation drop, same source depth), hence fluxes are antisymmetric
  // between neighbors: the update preserves positivity (h >= 0) and is
  // strictly mass-conservative, which makes the explicit scheme stable even
  // when dt is locally too large (no negative-depth clamping artifacts).
  float cap = 0.25f / dt;

  float out = min(Mp * fmax(0.f, Hc - Hxp), h * cap) +
              min(Mm * fmax(0.f, Hc - Hxm), h * cap) +
              min(Myp * fmax(0.f, Hc - Hyp), h * cap) +
              min(Mym * fmax(0.f, Hc - Hym), h * cap);

  float in = min(Mp * fmax(0.f, Hxp - Hc), hxp * cap) +
             min(Mm * fmax(0.f, Hxm - Hc), hxm * cap) +
             min(Myp * fmax(0.f, Hyp - Hc), hyp * cap) +
             min(Mym * fmax(0.f, Hym - Hc), hym * cap);

  // explicit update (fmax only guards against round-off)
  TSET(h_out, i, j, fmax(0.f, h + dt * (in - out)));
}
)""
