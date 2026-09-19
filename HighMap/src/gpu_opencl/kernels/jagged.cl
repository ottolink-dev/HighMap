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
                   const float          gamma,
                   const float          shape_gamma,
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

  // pass 2: compute distance to cell edge
  float w = 1.f;
  if (shape_gamma > 0.f)
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
    w = pow_float(tau, shape_gamma);
  }

  const sampler_t sampler_unnorm = CLK_NORMALIZED_COORDS_FALSE |
                                   CLK_ADDRESS_CLAMP_TO_EDGE |
                                   CLK_FILTER_NEAREST;

  float val_curr = read_imagef(in, sampler_unnorm, g).x;

  float effective_exponent = 1.f + (gamma - 1.f) * w;
  float val_out = val_curr > 0.f ? pow_float(val_curr, effective_exponent)
                                 : val_curr;

  if (has_mask > 0)
  {
    float t = mask[index];
    val_out = lerp(val_curr, val_out, t);
  }

  write_imagef(out, g, val_out);
}
)""
