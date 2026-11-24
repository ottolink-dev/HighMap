R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel generate_riverbed(global float   *sdf,
                              global float   *dz,
                              global float   *noise_x,
                              global float   *noise_y,
                              global float   *noise_r,
                              constant float *xp,
                              constant float *yp,
                              const int       nx,
                              const int       ny,
                              const int       np,
                              const float     depth_start,
                              const float     depth_end,
                              const float     slope_start,
                              const float     slope_end,
                              const float     shape_exponent_start,
                              const float     shape_exponent_end,
                              const float     k_smoothing,
                              const int       has_noise_x,
                              const int       has_noise_y,
                              const int       has_noise_r,
                              const float4    bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;
  float dr = has_noise_r > 0 ? noise_r[index] : 0.f;

  // (x, y) position
  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, dx, dy, bbox);

  // --- riverbed

  // relative position
  float2 e = (float2)(xp[np - 1] - xp[0], yp[np - 1] - yp[0]);
  float  dot_ee = dot(e, e);

  float2 w = (float2)(pos.x - xp[0], pos.y - yp[0]);
  float  t = clamp(dot(w, e) / dot_ee, 0.f, 1.f);

  float alpha = t * shape_exponent_start + (1.f - t) * shape_exponent_end;

  // radial function
  float r = sdf[index] + dr;
  r *= t * slope_start + (1.f - t) * slope_end;

  float dp = pow_float((pow_float(r, alpha) + 1.f), 1.f / alpha) - 1.f;
  dp = min_smooth(1.f, dp, k_smoothing);

  // adjust depth with distance
  float depth = t * depth_start + (1.f - t) * depth_end;
  dz[index] = (dp - 1.f) * depth;
}
)""
