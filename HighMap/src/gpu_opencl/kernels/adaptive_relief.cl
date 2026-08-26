R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel adaptive_relief(read_only image2d_t  in,
                            write_only image2d_t out,
                            const int            nx,
                            const int            ny,
                            const float          strength,
                            const float          clamp_ratio)
{
  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  float c = TGET(in, g.x, g.y);
  float n = TGET(in, g.x, g.y - 1);
  float s = TGET(in, g.x, g.y + 1);
  float w = TGET(in, g.x - 1, g.y);
  float e = TGET(in, g.x + 1, g.y);
  float nw = TGET(in, g.x - 1, g.y - 1);
  float ne = TGET(in, g.x + 1, g.y - 1);
  float sw = TGET(in, g.x - 1, g.y + 1);
  float se = TGET(in, g.x + 1, g.y + 1);

  // Local min / max
  float min_z = fmin(
      c,
      fmin(fmin(fmin(n, s), fmin(w, e)), fmin(fmin(nw, ne), fmin(sw, se))));
  float max_z = fmax(
      c,
      fmax(fmax(fmax(n, s), fmax(w, e)), fmax(fmax(nw, ne), fmax(sw, se))));

  // Discrete Laplacian (high-pass filter)
  float laplacian = c - (0.5f * (n + s + w + e) + 0.25f * (nw + ne + sw + se)) /
                            3.0f;

  float contrast = max_z - min_z;

  // Edge / slope magnitude
  float gx = (e - w) * 0.5f;
  float gy = (s - n) * 0.5f;
  float grad_sq = gx * gx + gy * gy;

  // Edge dampening factor
  float edge_factor = 1.0f /
                      (1.0f + 4.0f * grad_sq / (contrast * contrast + 1e-7f));

  float w_adaptive = strength * edge_factor;

  float val = c + w_adaptive * laplacian;

  if (clamp_ratio > 0.0f)
  {
    float clamped_val = clamp(val, min_z, max_z);
    val = mix(val, clamped_val, clamp_ratio);
  }

  TSET(out, g.x, g.y, (float4)(val, 0.0f, 0.0f, 1.0f));
}
)""
