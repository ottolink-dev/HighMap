R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
__constant const int   thermal_schott_di[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
__constant const int   thermal_schott_dj[8] = {0, 1, -1, 0, -1, 1, -1, 1};
__constant const float thermal_schott_c[8] =
    {1.f, 1.f, 1.f, 1.f, 1.41421356f, 1.41421356f, 1.41421356f, 1.41421356f};

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
    float dz = (val - z_in[linear_index(g.x + thermal_schott_di[k],
                                        g.y + thermal_schott_dj[k],
                                        nx)]) /
               thermal_schott_c[k];

    if (dz > talus_val)
      down++;
    else if (dz < -talus_val)
      up++;
  }

  z_out[index] = val + intensity * talus_val * (float)(up - down);
}
)""
