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
                            const float4    bbox,
                            const int       is_closed)
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
  float s = 1.f;

  for (int k = 0; k < np - 1; ++k)
  {
    int    p = k + 1;
    float2 p_k = (float2)(xp[k], yp[k]);
    float2 p_next = (float2)(xp[p], yp[p]);

    float2 e = p_next - p_k;
    float2 w = pos - p_k;
    float  t = clamp(dot(w, e) / dot(e, e), 0.f, 1.f);
    float2 b = w - e * t;
    d = min(d, dot(b, b));

    if (is_closed > 0)
    {
      int c0 = pos.y >= p_k.y ? 1 : 0;
      int c1 = pos.y < p_next.y ? 1 : 0;
      int c2 = (e.x * w.y > e.y * w.x) ? 1 : 0;

      if ((c0 && c1 && c2) || (!c0 && !c1 && !c2)) s *= -1.f;
    }
  }

  sdf2[index] = s * sqrt(d);
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
                                   const float4    bbox,
                                   const int       is_closed)
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
  int   crossings = 0;

  for (int i = 0; i < np - 1; i += 3)
  {
    int j = i + 1;
    int k = i + 2;

    float2 p0 = (float2)(xp[i], yp[i]);
    float2 p1 = (float2)(xp[j], yp[j]);
    float2 p2 = (float2)(xp[k], yp[k]);

    float2 a = p1 - p0;
    float2 b = p0 - 2.f * p1 + p2;
    float2 c = 2.f * a;
    float2 d = p0 - pos;

    float b_len_sq = dot(b, b);
    float dist = 0.f;

    if (b_len_sq < 1e-10f)
    {
      // Colinear degenerate quadratic bezier: distance to segment [p0, p2]
      float2 e = p2 - p0;
      float2 w = pos - p0;
      float  lsq = dot(e, e);
      float  t = (lsq > 0.f) ? clamp(dot(w, e) / lsq, 0.f, 1.f) : 0.f;
      float2 diff = w - e * t;
      dist = dot(diff, diff);
    }
    else
    {
      float kk = 1.f / b_len_sq;
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
        float  zz = sqrt(max(0.f, -p));
        float  denom = p * zz * 2.f;
        float  arg = (fabs(denom) > 1e-12f) ? clamp(q / denom, -1.f, 1.f) : 0.f;
        float  v = acos(arg) / 3.f;
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
    }

    dmin = min(dmin, dist);

    // Closed path ray-crossing parity test
    if (is_closed > 0)
    {
      float a_y = p0.y - 2.f * p1.y + p2.y;
      float b_y = 2.f * (p1.y - p0.y);
      float c_y = p0.y - pos.y;

      if (fabs(a_y) < 1e-9f)
      {
        if (fabs(b_y) > 1e-9f)
        {
          float t = -c_y / b_y;
          if (t >= 0.f && t < 1.f)
          {
            float bx = (1.f - t) * (1.f - t) * p0.x +
                       2.f * t * (1.f - t) * p1.x + t * t * p2.x;
            if (bx > pos.x) crossings++;
          }
        }
      }
      else
      {
        float discr = b_y * b_y - 4.f * a_y * c_y;
        if (discr >= 0.f)
        {
          float s_discr = sqrt(discr);
          float t1 = (-b_y - s_discr) / (2.f * a_y);
          float t2 = (-b_y + s_discr) / (2.f * a_y);

          if (t1 >= 0.f && t1 < 1.f)
          {
            float bx1 = (1.f - t1) * (1.f - t1) * p0.x +
                        2.f * t1 * (1.f - t1) * p1.x + t1 * t1 * p2.x;
            if (bx1 > pos.x) crossings++;
          }
          if (t2 >= 0.f && t2 < 1.f)
          {
            float bx2 = (1.f - t2) * (1.f - t2) * p0.x +
                        2.f * t2 * (1.f - t2) * p1.x + t2 * t2 * p2.x;
            if (bx2 > pos.x) crossings++;
          }
        }
      }
    }
  }

  float s = (is_closed > 0 && (crossings % 2 == 1)) ? -1.f : 1.f;
  sdf2[index] = s * sqrt(dmin);
}
)""
