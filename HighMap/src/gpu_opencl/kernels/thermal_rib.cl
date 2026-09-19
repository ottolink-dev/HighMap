R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

__constant const int   thermal_rib_di[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
__constant const int   thermal_rib_dj[8] = {0, 1, -1, 0, -1, 1, -1, 1};
__constant const float thermal_rib_c[8] =
    {1.f, 1.f, 1.f, 1.f, 1.41421356f, 1.41421356f, 1.41421356f, 1.41421356f};

void kernel thermal_rib(global const float *z_in,
                        global float       *z_out,
                        const int           nx,
                        const int           ny)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  int  index = linear_index(g.x, g.y, nx);

  if (g.x >= nx || g.y >= ny) return;
  if (apply_boundaries_io(z_in, z_out, g.x, g.y, nx, ny)) return;

  // --- thermal erosion

  float delta_min = FLT_MAX;
  float val = z_in[index];

#pragma unroll
  for (int k = 0; k < 8; k++)
  {
    float dz = fabs(val - z_in[linear_index(g.x + thermal_rib_di[k],
                                            g.y + thermal_rib_dj[k],
                                            nx)]) /
               thermal_rib_c[k];
    delta_min = min(delta_min, dz);
  }

  z_out[index] = val - delta_min;
}
)""
