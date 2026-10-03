R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel voronoi_shrink(read_only image2d_t  in,
                           write_only image2d_t out,
                           global const float  *mask,
                           global const float  *noise_x,
                           global const float  *noise_y,
                           const int            nx,
                           const int            ny,
                           const float          kx,
                           const float          ky,
                           const uint           seed,
                           const float2         jitter,
                           const float          shrink_factor,
                           const float          fill_value,
                           const float          angle,
                           const int            has_mask,
                           const int            has_noise_x,
                           const int            has_noise_y,
                           const float4         bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  float alpha = angle / 180.f * 3.14159265f;
  float ca = cos(alpha);
  float sa = sin(alpha);

  if (angle != 0.f)
  {
    pos = (float2)(pos.x * ca - pos.y * sa, pos.x * sa + pos.y * ca);
  }

  float2 p = floor(pos);
  float2 pi;
  float2 f = fract(pos, &pi);

  float2 mb = (float2)(0.f, 0.f);
  float2 mr = (float2)(0.f, 0.f);
  float  res = 8.f;

  // pass 1: find closest Voronoi cell center
  for (int dy_i = -2; dy_i <= 2; dy_i++)
    for (int dx_i = -2; dx_i <= 2; dx_i++)
    {
      float2 b = (float2)(dx_i, dy_i);
      float2 r = b - f + jitter * hash22f(p + b, fseed);
      float  d = dot(r, r);

      if (d < res)
      {
        res = d;
        mr = r;
        mb = b;
      }
    }

  // cell center in continuous grid coordinates
  float2 center = p + mb + jitter * hash22f(p + mb, fseed);

  // linear shrink towards center
  float  s = max(shrink_factor, 1e-5f);
  float2 pos_src = center + (pos - center) / s;

  // verify if pos_src belongs to the same Voronoi cell
  float2 p_src = floor(pos_src);
  float2 pi_src;
  float2 f_src = fract(pos_src, &pi_src);

  float2 mb_src = (float2)(0.f, 0.f);
  float  res_src = 8.f;

  for (int dy_i = -2; dy_i <= 2; dy_i++)
    for (int dx_i = -2; dx_i <= 2; dx_i++)
    {
      float2 b = (float2)(dx_i, dy_i);
      float2 r = b - f_src + jitter * hash22f(p_src + b, fseed);
      float  d = dot(r, r);

      if (d < res_src)
      {
        res_src = d;
        mb_src = b;
      }
    }

  float2 center_src = p_src + mb_src + jitter * hash22f(p_src + mb_src, fseed);

  bool in_same_cell = (dot(center - center_src, center - center_src) < 1e-5f);

  const sampler_t sampler_curr = CLK_NORMALIZED_COORDS_FALSE |
                                 CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const sampler_t sampler_interp = CLK_NORMALIZED_COORDS_TRUE |
                                   CLK_ADDRESS_CLAMP_TO_EDGE |
                                   CLK_FILTER_LINEAR;

  float val_curr = read_imagef(in, sampler_curr, g).x;
  float val_out;

  if (in_same_cell)
  {
    // unrotate if needed
    float2 pos_unrot = pos_src;
    if (angle != 0.f)
    {
      pos_unrot = (float2)(pos_src.x * ca + pos_src.y * sa,
                           -pos_src.x * sa + pos_src.y * ca);
    }

    // map from Voronoi domain back to normalized grid coordinates [0, 1]
    float u = (pos_unrot.x / kx - bbox.x) / (bbox.y - bbox.x);
    float v = (pos_unrot.y / ky - bbox.z) / (bbox.w - bbox.z);

    // in g_to_xy: x = (float)g.x / (float)nx, so g.x / nx = u.
    // In OpenCL normalized coordinates, pixel g.x center is (g.x + 0.5) / nx =
    // u + 0.5 / nx.
    float2 coord = (float2)(u + 0.5f / (float)nx, v + 0.5f / (float)ny);

    val_out = read_imagef(in, sampler_interp, coord).x;
  }
  else
  {
    val_out = fill_value;
  }

  if (has_mask > 0)
  {
    float t = mask[index];
    val_out = lerp(val_curr, val_out, t);
  }

  write_imagef(out, g, val_out);
}

void kernel voronoi_shrink_points(read_only image2d_t  in,
                                  write_only image2d_t out,
                                  global const float  *mask,
                                  global const float  *noise_x,
                                  global const float  *noise_y,
                                  constant float      *xp,
                                  constant float      *yp,
                                  const int            nx,
                                  const int            ny,
                                  const int            np,
                                  const float          shrink_factor,
                                  const float          fill_value,
                                  const int            has_mask,
                                  const int            has_noise_x,
                                  const int            has_noise_y,
                                  const float4         bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, dx, dy, bbox);

  // find closest point in cloud
  float dmin = FLT_MAX;
  int   idx_closest = 0;

  for (int i = 0; i < np; ++i)
  {
    float2 pt = (float2)(xp[i], yp[i]);
    float2 diff = pos - pt;
    float  d = dot(diff, diff);
    if (d < dmin)
    {
      dmin = d;
      idx_closest = i;
    }
  }

  float2 center = (float2)(xp[idx_closest], yp[idx_closest]);

  // linear shrink towards center
  float  s = max(shrink_factor, 1e-5f);
  float2 pos_src = center + (pos - center) / s;

  // verify if pos_src belongs to the same Voronoi cell
  float dmin_src = FLT_MAX;
  int   idx_closest_src = 0;

  for (int i = 0; i < np; ++i)
  {
    float2 pt = (float2)(xp[i], yp[i]);
    float2 diff = pos_src - pt;
    float  d = dot(diff, diff);
    if (d < dmin_src)
    {
      dmin_src = d;
      idx_closest_src = i;
    }
  }

  bool in_same_cell = (idx_closest == idx_closest_src);

  const sampler_t sampler_curr = CLK_NORMALIZED_COORDS_FALSE |
                                 CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const sampler_t sampler_interp = CLK_NORMALIZED_COORDS_TRUE |
                                   CLK_ADDRESS_CLAMP_TO_EDGE |
                                   CLK_FILTER_LINEAR;

  float val_curr = read_imagef(in, sampler_curr, g).x;
  float val_out;

  if (in_same_cell)
  {
    // map from domain coordinates back to normalized grid coordinates [0, 1]
    float u = (pos_src.x - bbox.x) / (bbox.y - bbox.x);
    float v = (pos_src.y - bbox.z) / (bbox.w - bbox.z);

    float2 coord = (float2)(u + 0.5f / (float)nx, v + 0.5f / (float)ny);

    val_out = read_imagef(in, sampler_interp, coord).x;
  }
  else
  {
    val_out = fill_value;
  }

  if (has_mask > 0)
  {
    float t = mask[index];
    val_out = lerp(val_curr, val_out, t);
  }

  write_imagef(out, g, val_out);
}
)""
