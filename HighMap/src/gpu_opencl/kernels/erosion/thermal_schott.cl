R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel thermal_schott(global const float *z_in,
                           global float       *z_out,
                           const global float *talus,
                           const int           nx,
                           const int           ny,
                           const float         intensity)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  int  index = linear_index(g.x, g.y, nx);

  if (g.x >= nx || g.y >= ny) return;
  if (apply_boundaries_io(z_in, z_out, g.x, g.y, nx, ny)) return;

  // --- thermal erosion

  float val = z_in[index];
  float talus_val = talus[index];

  int up = 0;
  int down = 0;

#pragma unroll
  for (int k = 0; k < 8; k++)
  {
    float dz = (val - z_in[linear_index(g.x + d8_di[k], g.y + d8_dj[k], nx)]) /
               d8_dist[k];

    if (dz > talus_val)
      down++;
    else if (dz < -talus_val)
      up++;
  }

  z_out[index] = val + intensity * talus_val * (float)(up - down);
}
)""
