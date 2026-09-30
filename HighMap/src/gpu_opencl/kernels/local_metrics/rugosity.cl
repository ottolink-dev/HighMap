R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel rugosity_post(global float *z_skw,
                          global float *z_std,
                          const int     nx,
                          const int     ny,
                          const float   tol,
                          const int     convex)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int   index = linear_index(g.x, g.y, nx);
  float skw = z_skw[index];
  float std = z_std[index];

  if (std > tol) skw /= pow_float(std, 1.5f);

  if (convex > 0)
    skw = max(0.f, skw);
  else
    skw = min(0.f, skw);

  z_skw[index] = skw;
}
)""
