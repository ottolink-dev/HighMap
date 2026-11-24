R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel laplace(global float *array,
                    const int     nx,
                    const int     ny,
                    const float   sigma)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int im1 = g.x == 0 ? 1 : g.x - 1;
  int ip1 = g.x == nx - 1 ? nx - 2 : g.x + 1;
  int jm1 = g.y == 0 ? 1 : g.y - 1;
  int jp1 = g.y == ny - 1 ? ny - 2 : g.y + 1;

  // compute laplacian
  float delta = -4.f * array[linear_index(g.x, g.y, nx)] +
                array[linear_index(ip1, g.y, nx)] +
                array[linear_index(im1, g.y, nx)] +
                array[linear_index(g.x, jm1, nx)] +
                array[linear_index(g.x, jp1, nx)];

  array[linear_index(g.x, g.y, nx)] += sigma * delta;
}
)""
