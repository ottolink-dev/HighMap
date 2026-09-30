R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel project_talus_along_direction(global float *array,
                                          global float *talus,
                                          global float *out,
                                          const int     nx,
                                          const int     ny,
                                          const int     di,
                                          const int     dj,
                                          const float   vmin)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int   index = linear_index(g.x, g.y, nx);
  float val = array[index];

  if (val < vmin)
  {
    out[index] = val;
    return;
  }

  // full domain diagonal in worst case scenario
  const int nsteps = (int)(1.414f * max(nx, ny));

  for (int k = 0; k < nsteps; ++k)
  {
    int   i = g.x + (k + 1) * di;
    int   j = g.y + (k + 1) * dj;
    float dr = length((float2)(i - g.x, j - g.y));

    if (i < 0 || i >= nx || j < 0 || j >= ny) break;

    int   idx = linear_index(i, j, nx);
    float v = val - dr * talus[idx];
    float vref = out[idx];

    if (v < vref) break;

    atomic_max_float(&out[idx], v);
  }
}
)""
