R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel
harmonic_interpolation_red_black(global float       *u,
                                 const global float *mask_fixed_values,
                                 const int           nx,
                                 const int           ny,
                                 const int           color,
                                 const float         omega)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x <= 0 || g.x >= nx - 1 || g.y <= 0 || g.y >= ny - 1) return;

  // Process only cells matching the active color (0 = Red, 1 = Black)
  if (((g.x + g.y) & 1) != color) return;

  int idx = linear_index(g.x, g.y, nx);
  if (mask_fixed_values[idx] > 0.f) return;

  float u_old = u[idx];
  float u_avg = 0.25f * (u[linear_index(g.x - 1, g.y, nx)] +
                         u[linear_index(g.x + 1, g.y, nx)] +
                         u[linear_index(g.x, g.y - 1, nx)] +
                         u[linear_index(g.x, g.y + 1, nx)]);

  u[idx] = u_old + omega * (u_avg - u_old);
}

void kernel
harmonic_interpolation_red_black_diff(global float          *u,
                                      const global float    *mask_fixed_values,
                                      volatile global float *g_max_diff,
                                      const int              nx,
                                      const int              ny,
                                      const int              color,
                                      const float            omega,
                                      const int              track_diff)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x <= 0 || g.x >= nx - 1 || g.y <= 0 || g.y >= ny - 1) return;

  if (((g.x + g.y) & 1) != color) return;

  int idx = linear_index(g.x, g.y, nx);
  if (mask_fixed_values[idx] > 0.f) return;

  float u_old = u[idx];
  float u_avg = 0.25f * (u[linear_index(g.x - 1, g.y, nx)] +
                         u[linear_index(g.x + 1, g.y, nx)] +
                         u[linear_index(g.x, g.y - 1, nx)] +
                         u[linear_index(g.x, g.y + 1, nx)]);

  float diff = omega * (u_avg - u_old);
  u[idx] = u_old + diff;

  if (track_diff > 0)
  {
    atomic_max_float(g_max_diff, fabs(diff));
  }
}

void kernel harmonic_interpolation_red_black_diff_aniso(
    global float          *u,
    const global float    *mask_fixed_values,
    const global float    *dx,
    const global float    *dy,
    volatile global float *g_max_diff,
    const int              nx,
    const int              ny,
    const int              color,
    const float            omega,
    const int              track_diff)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x <= 0 || g.x >= nx - 1 || g.y <= 0 || g.y >= ny - 1) return;

  if (((g.x + g.y) & 1) != color) return;

  int idx = linear_index(g.x, g.y, nx);
  if (mask_fixed_values[idx] > 0.f) return;

  int idx_w = linear_index(g.x - 1, g.y, nx);
  int idx_e = linear_index(g.x + 1, g.y, nx);
  int idx_s = linear_index(g.x, g.y - 1, nx);
  int idx_n = linear_index(g.x, g.y + 1, nx);

  float dx_c = dx[idx];
  float dy_c = dy[idx];

  float k_w = 0.5f * (dx[idx_w] + dx_c);
  float k_e = 0.5f * (dx[idx_e] + dx_c);
  float k_s = 0.5f * (dy[idx_s] + dy_c);
  float k_n = 0.5f * (dy[idx_n] + dy_c);

  float c_tot = k_w + k_e + k_s + k_n;
  if (c_tot < 1e-8f) return;

  float u_old = u[idx];
  float u_star = (k_w * u[idx_w] + k_e * u[idx_e] + k_s * u[idx_s] +
                  k_n * u[idx_n]) /
                 c_tot;

  float diff = omega * (u_star - u_old);
  u[idx] = u_old + diff;

  if (track_diff > 0)
  {
    atomic_max_float(g_max_diff, fabs(diff));
  }
}
)""
