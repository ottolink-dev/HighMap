R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel phase_averaging(read_only image2d_t  fr,
                            read_only image2d_t  fi,
                            write_only image2d_t fr_out,
                            write_only image2d_t fi_out,
                            int                  nx,
                            int                  ny,
                            int                  ir)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int i = g.x;
  int j = g.y;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float2 sum = 0.f;
  float2 wsum = 0.f;
  float  dmax = (float)ir;

  for (int p = -ir; p <= ir; ++p)
    for (int q = -ir; q <= ir; ++q)
    {
      float r = length((float2)(p, q)) / dmax;

      if (r > 1.f) continue;

      float w = 1.f - r * r * (3.f - 2.f * r);
      sum += w * (float2)(TGET(fr, i + p, j + q), TGET(fi, i + p, j + q));
      // count++;
    }

  // val /= count;
  // val /= length(val);

  TSET(fr_out, i, j, sum.x);
  TSET(fi_out, i, j, sum.y);
}
)""
