R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

/* Port of erosiv/geotransport (MIT, McDonald & Cordonnier): stochastic
 * Monte-Carlo estimator for steady-state transport / flow accumulation. */

float facc_st_bilinear(global const float *a, float x, float y, int nx, int ny)
{
  int i = (int)x;
  int j = (int)y;
  if (i > nx - 2) i = nx - 2;
  if (j > ny - 2) j = ny - 2;
  float u = x - (float)i;
  float v = y - (float)j;
  return (1.f - u) * (1.f - v) * a[linear_index(i, j, nx)] +
         u * (1.f - v) * a[linear_index(i + 1, j, nx)] +
         (1.f - u) * v * a[linear_index(i, j + 1, nx)] +
         u * v * a[linear_index(i + 1, j + 1, nx)];
}

float facc_st_stepsize(float px, float py, float dx, float dy)
{
  const float tmax = 1.41421356f;

  float x_neg = floor(px);
  float y_neg = floor(py);

  float tx = min(max((x_neg - px) / dx, (1.f + x_neg - px) / dx), tmax);
  float ty = min(max((y_neg - py) / dy, (1.f + y_neg - py) / dy), tmax);

  return 0.5f * (tx + ty);
}

void kernel flow_accum_stochastic_solve(global const float *vx,
                                        global const float *vy,
                                        global const float *source,
                                        global const float *decay,
                                        global float       *flux,
                                        const int           has_source,
                                        const int           has_decay,
                                        const int           nx,
                                        const int           ny,
                                        const int           n_samples,
                                        const uint          seed)
{
  int n = get_global_id(0);
  if (n >= n_samples) return;

  const float eps = 1e-16f;
  const int   max_step = nx + ny;
  const float p_inv = (float)nx * (float)ny;

  float px = (float)(hash21u((uint)n, seed) >> 8) / 16777216.f * (float)nx;
  float py = (float)(hash21u((uint)n, seed + 0x9e3779b9u) >> 8) / 16777216.f *
             (float)ny;

  int i0 = (int)px;
  int j0 = (int)py;
  if (i0 >= nx || j0 >= ny) return;

  float s = (has_source > 0 ? source[linear_index(i0, j0, nx)] : 1.f) * p_inv;
  if (fabs(s) < eps) return;

  float att = 1.f;
  int   ind = linear_index(i0, j0, nx);

  for (int step = 0; step < max_step; ++step)
  {
    int ci = (int)px;
    int cj = (int)py;
    if (ci < 0 || ci >= nx || cj < 0 || cj >= ny) break;

    int cind = linear_index(ci, cj, nx);
    if (cind != ind)
    {
      ind = cind;
      atomic_add_float(&flux[cind], s * att);
    }

    float v_x = facc_st_bilinear(vx, px, py, nx, ny);
    float v_y = facc_st_bilinear(vy, px, py, nx, ny);
    float v_len = hypot(v_x, v_y);
    if (v_len < eps) break;

    float ds = facc_st_stepsize(px, py, v_x / v_len, v_y / v_len);
    px += ds * v_x / v_len;
    py += ds * v_y / v_len;

    float dec = has_decay > 0 ? decay[linear_index(ci, cj, nx)] : 0.f;
    att *= exp(-ds * 1.41421356f / v_len * dec);
    if (fabs(att) < eps) break;
  }
}

void kernel flow_accum_stochastic_normalize(global float       *flux,
                                            global const float *vx,
                                            global const float *vy,
                                            global const float *source,
                                            const int           has_source,
                                            const int           nx,
                                            const int           ny,
                                            const int           n_samples)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  if (g.x >= nx || g.y >= ny) return;

  int   idx = linear_index(g.x, g.y, nx);
  float src = has_source > 0 ? source[idx] : 1.f;
  float norm = fabs(vx[idx]) + fabs(vy[idx]);
  norm = max(norm, 1e-9f);

  flux[idx] = (src + flux[idx] / (float)n_samples) / norm;
}
)""
