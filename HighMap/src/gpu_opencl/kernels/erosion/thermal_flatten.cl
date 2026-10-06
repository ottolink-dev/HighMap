R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel thermal_flatten(global const float *z_in,
                            global float       *z_out,
                            const global float *talus,
                            const int           nx,
                            const int           ny,
                            const float         sigma_inf,
                            const float         sigma_sup)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  int  index = linear_index(g.x, g.y, nx);

  if (g.x >= nx || g.y >= ny) return;
  if (apply_boundaries_buffer_io(z_in, z_out, g.x, g.y, nx, ny, 2)) return;

  const int   rotation_offset = (g.x + g.y) & 7;
  const float zc = z_in[index];

  const float slope_min = 0.f;

  // --- Remove material from current cell

  float delta = 0.f;
  float dmax = 0.f;
  int   ka = -1;

#pragma unroll
  for (int n = 0; n < 8; ++n)
  {
    const int   k = (n + rotation_offset) & 7;
    const int   ni = g.x + d8_di[k];
    const int   nj = g.y + d8_dj[k];
    const float dz = (zc - z_in[linear_index(ni, nj, nx)]) / d8_dist[k];

    if (dz > dmax)
    {
      dmax = dz;
      ka = k;
    }
  }

  float excess = max(dmax - slope_min, 0.f);

  if (dmax < talus[index]) delta -= sigma_inf * excess;
  if (dmax > talus[index]) delta -= sigma_sup * excess;

  // --- Remove material from neighbors

#pragma unroll
  for (int n = 0; n < 8; ++n)
  {
    const int   k = (n + rotation_offset) & 7;
    const int   pi = g.x - d8_di[k];
    const int   pj = g.y - d8_dj[k];
    const int   pidx = linear_index(pi, pj, nx);
    const float zp = z_in[pidx];

    float dmax_p = 0.f;
    int   ka_p = -1;

// find steepest neighbor for source cell
#pragma unroll
    for (int m = 0; m < 8; ++m)
    {
      const int   kk = (m + rotation_offset) & 7;
      const int   qi = pi + d8_di[kk];
      const int   qj = pj + d8_dj[kk];
      const float dz = (zp - z_in[linear_index(qi, qj, nx)]) / d8_dist[kk];

      if (dz > dmax_p)
      {
        dmax_p = dz;
        ka_p = kk;
      }
    }

    // if current cell is target receiver
    if (ka_p == k && dmax_p > 0.0f && dmax_p < talus[pidx])
      delta += 0.5f * dmax_p;
  }

  z_out[index] = zc + delta;
}
)""
