R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
__constant const int   thermal_scree_di[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
__constant const int   thermal_scree_dj[8] = {0, 1, -1, 0, -1, 1, -1, 1};
__constant const float thermal_scree_c[8] =
    {1.f, 1.f, 1.f, 1.f, 1.41421356f, 1.41421356f, 1.41421356f, 1.41421356f};

void kernel thermal_scree(global const float *z_in,
                          global float       *z_out,
                          const global float *talus,
                          const global float *zmax,
                          const int           nx,
                          const int           ny)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  int  index = linear_index(g.x, g.y, nx);

  if (g.x >= nx || g.y >= ny) return;
  if (apply_boundaries_io(z_in, z_out, g.x, g.y, nx, ny)) return;

  // --- thermal erosion

  float val = z_in[index];
  float sum = 0.f;
  float slope_max = 0.f;

#pragma unroll
  for (int k = 0; k < 8; k++)
  {
    float dz = (val - z_in[linear_index(g.x + thermal_scree_di[k],
                                        g.y + thermal_scree_dj[k],
                                        nx)]) /
               thermal_scree_c[k];

    if (dz < 0.f) sum += dz;

    slope_max = max(slope_max, fabs(dz));
  }

  float t = talus[index];
  float amp = clamp(1.f - t / slope_max, 0.f, 1.f);
  amp = smoothstep3(amp);

  amp *= almost_unit_identity(clamp(1.f - val / zmax[index], 0.f, 1.f));

  z_out[index] = val + 0.25f * (t - 0.5f * sum) * amp;
}
)""
