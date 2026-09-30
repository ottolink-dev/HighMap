R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel local_skewness(read_only image2d_t  img_in,
                           write_only image2d_t img_out,
                           const int            nx,
                           const int            ny,
                           const int            ir)
{
  const int2 g = {get_global_id(0), get_global_id(1)};
  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float s = 0.f;
  float ss = 0.f;
  float sss = 0.f;
  int   count = 0;

  for (int i = g.x - ir; i <= g.x + ir; ++i)
    for (int j = g.y - ir; j <= g.y + ir; ++j)
      if (is_within_radius_and_inside(i, j, g.x, g.y, ir, nx, ny))
      {
        float v = read_imagef(img_in, sampler, (int2)(i, j)).x;
        s += v;
        ss += v * v;
        sss += v * v * v;
        count++;
      }

  if (count < 2)
  {
    write_imagef(img_out, g, 0.f);
    return;
  }

  float mean = s / (float)count;

  float var = max(0.f, (count * ss - s * s)) /
              ((float)count * (float)(count - 1));

  float std = sqrt(var);
  float m3 = (sss - 3.f * mean * ss + 2.f * mean * mean * s) / (float)count;
  float skew = (std > 1e-12f) ? (m3 / (std * std * std)) : 0.f;

  write_imagef(img_out, g, skew);
}
)""
