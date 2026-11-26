R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel poisson_compute_laplacian(const global float *src,
                                      global float       *laplacian,
                                      const int           nx,
                                      const int           ny)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int im1 = g.x == 0 ? 1 : g.x - 1;
  int ip1 = g.x == nx - 1 ? nx - 2 : g.x + 1;
  int jm1 = g.y == 0 ? 1 : g.y - 1;
  int jp1 = g.y == ny - 1 ? ny - 2 : g.y + 1;

  float delta = -4.f * src[linear_index(g.x, g.y, nx)] +
                src[linear_index(ip1, g.y, nx)] +
                src[linear_index(im1, g.y, nx)] +
                src[linear_index(g.x, jm1, nx)] +
                src[linear_index(g.x, jp1, nx)];

  laplacian[linear_index(g.x, g.y, nx)] = delta;
}

void kernel blend_poisson_red_black(global float       *array1,
                                    const global float *delta2,
                                    const global float *mask,
                                    const int           nx,
                                    const int           ny,
                                    const int           color,
                                    const float         omega,
                                    const int           has_mask)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;
  if (((g.x + g.y) & 1) != color) return;

  int im1 = g.x == 0 ? 1 : g.x - 1;
  int ip1 = g.x == nx - 1 ? nx - 2 : g.x + 1;
  int jm1 = g.y == 0 ? 1 : g.y - 1;
  int jp1 = g.y == ny - 1 ? ny - 2 : g.y + 1;

  int idx = linear_index(g.x, g.y, nx);

  float u_curr = array1[idx];
  float sum_nb = array1[linear_index(ip1, g.y, nx)] +
                 array1[linear_index(im1, g.y, nx)] +
                 array1[linear_index(g.x, jm1, nx)] +
                 array1[linear_index(g.x, jp1, nx)];

  float rhs = delta2[idx];
  float u_star = 0.25f * (sum_nb - rhs);

  float amp = has_mask > 0 ? mask[idx] : 1.f;

  array1[idx] = u_curr + omega * amp * (u_star - u_curr);
}
)""
