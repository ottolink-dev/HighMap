R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel jagged(read_only image2d_t  in,
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
                   const float          factor,
                   const float          shape_factor,
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

  float2 p = floor(pos);
  float2 pi;
  float2 f = fract(pos, &pi);

  float2 mb = (float2)(0.f, 0.f);
  float2 mr = (float2)(0.f, 0.f);
  float  res = 8.f;

  // pass 1: find closest Voronoi cell center
  for (int dy_i = -1; dy_i <= 1; dy_i++)
    for (int dx_i = -1; dx_i <= 1; dx_i++)
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

  float2 best_feature_point = pos + mr;

  // pass 2: compute distance to cell edge if shape factor > 0
  float w = 1.f;
  if (shape_factor > 0.f)
  {
    float d_edge = 8.f;
    for (int dy_i = -2; dy_i <= 2; dy_i++)
      for (int dx_i = -2; dx_i <= 2; dx_i++)
      {
        float2 b = mb + (float2)(dx_i, dy_i);
        float2 r = b - f + jitter * hash22f(p + b, fseed);
        if (dot(mr - r, mr - r) > 1e-5f)
        {
          float d = dot(0.5f * (mr + r), normalize(r - mr));
          d_edge = min(d_edge, d);
        }
      }
    d_edge = max(0.f, d_edge);
    float d_center = length(mr);
    float tau = d_edge / (d_center + d_edge + 1e-6f);
    tau = smoothstep3(clamp(tau, 0.f, 1.f));
    w = pow_float(tau, shape_factor);
  }

  // map best_feature_point in pos-space back to normalized domain coordinates
  // [0, 1]
  float u_c = (best_feature_point.x / kx - bbox.x) / (bbox.y - bbox.x);
  float v_c = (best_feature_point.y / ky - bbox.z) / (bbox.w - bbox.z);

  const sampler_t sampler_norm = CLK_NORMALIZED_COORDS_TRUE |
                                 CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  const sampler_t sampler_unnorm = CLK_NORMALIZED_COORDS_FALSE |
                                   CLK_ADDRESS_CLAMP_TO_EDGE |
                                   CLK_FILTER_NEAREST;

  float val_voronoi = read_imagef(in, sampler_norm, (float2)(u_c, v_c)).x;
  float val_curr = read_imagef(in, sampler_unnorm, g).x;

  float val_jagged = val_voronoi + factor * (val_voronoi - val_curr);
  float val_out = lerp(val_curr, val_jagged, w);

  if (has_mask > 0)
  {
    float t = mask[index];
    val_out = lerp(val_curr, val_out, t);
  }

  write_imagef(out, g, val_out);
}
)""
