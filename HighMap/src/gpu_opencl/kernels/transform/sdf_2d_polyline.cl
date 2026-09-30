R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel sdf_2d_polyline(global float   *sdf2,
                            global float   *noise_x,
                            global float   *noise_y,
                            constant float *xp,
                            constant float *yp,
                            const int       nx,
                            const int       ny,
                            const int       np,
                            const int       has_noise_x,
                            const int       has_noise_y,
                            const float4    bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  // (x, y) position
  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, dx, dy, bbox);

  // sdf
  float d = FLT_MAX;

  for (int k = 0; k < np - 1; ++k)
  {
    int    p = k + 1;
    float2 e = (float2)(xp[p] - xp[k], yp[p] - yp[k]);
    float2 w = (float2)(pos.x - xp[k], pos.y - yp[k]);
    float  t = clamp(dot(w, e) / dot(e, e), 0.f, 1.f);
    float2 b = w - e * t;
    d = min(d, dot(b, b));
  }

  sdf2[index] = sqrt(d);
}

void kernel sdf_2d_polyline_bezier(global float   *sdf2,
                                   global float   *noise_x,
                                   global float   *noise_y,
                                   constant float *xp,
                                   constant float *yp,
                                   const int       nx,
                                   const int       ny,
                                   const int       np,
                                   const int       has_noise_x,
                                   const int       has_noise_y,
                                   const float4    bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  // (x, y) position
  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, dx, dy, bbox);

  // sdf
  float dmin = FLT_MAX;

  for (int i = 0; i < np - 1; i += 3)
  {
    int j = i + 1;
    int k = i + 2;

    float dist = 0.f;

    float2 a = (float2)(xp[j] - xp[i], yp[j] - yp[i]);
    float2 b = (float2)(xp[i] - 2.f * xp[j] + xp[k],
                        yp[i] - 2.f * yp[j] + yp[k]);
    float2 c = (float2)(2.f * a.x, 2.f * a.y);
    float2 d = (float2)(xp[i] - pos.x, yp[i] - pos.y);

    float kk = 1.f / dot(b, b);
    float kx = kk * dot(a, b);
    float ky = kk * (2.f * dot(a, a) + dot(d, b)) / 3.f;
    float kz = kk * dot(d, a);
    float p = ky - kx * kx;
    float p3 = p * p * p;
    float q = kx * (2.f * kx * kx - 3.f * ky) + kz;
    float h = q * q + 4.f * p3;

    if (h >= 0.f)
    {
      h = sqrt(h);
      float2 xx = (float2)(0.5f * (h - q), 0.5f * (-h - q));
      float2 uv = (float2)(copysign(1.f, xx.x) *
                               pow_float(fabs(xx.x), 1.f / 3.f),
                           copysign(1.f, xx.y) *
                               pow_float(fabs(xx.y), 1.f / 3.f));
      float  t = clamp(uv.x + uv.y - kx, 0.f, 1.f);

      float2 dd = (float2)(d.x + (c.x + b.x * t) * t,
                           d.y + (c.y + b.y * t) * t);

      dist = dot(dd, dd);
    }
    else
    {
      float  zz = sqrt(-p);
      float  v = acos(q / (p * zz * 2.f)) / 3.f;
      float  m = cos(v);
      float  n = sin(v) * 1.732050808f;
      float3 tt = (float3)(clamp((m + m) * zz - kx, 0.f, 1.f),
                           clamp((-n - m) * zz - kx, 0.f, 1.f),
                           clamp((n - m) * zz - kx, 0.f, 1.f));

      float2 dd1 = (float2)(d.x + (c.x + b.x * tt.x) * tt.x,
                            d.y + (c.y + b.y * tt.x) * tt.x);
      float2 dd2 = (float2)(d.x + (c.x + b.x * tt.y) * tt.y,
                            d.y + (c.y + b.y * tt.y) * tt.y);

      float dist1 = dot(dd1, dd1);
      float dist2 = dot(dd2, dd2);

      dist = min(dist1, dist2);
    }

    dmin = min(dmin, dist);
  }

  sdf2[index] = sqrt(dmin);
}
)""
