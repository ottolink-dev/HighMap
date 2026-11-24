R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel laplacian_fract(read_only image2d_t  img_in,
                            write_only image2d_t img_out,
                            const int            nx,
                            const int            ny,
                            const int            ir,
                            const float          s)
{
  // Adapted from terrain-descriptors
  // Original author: oargudo
  // License: MIT (see THIRD_PARTY_LICENSES.md)

  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float sum = 0.f;

  const float h_center = read_imagef(img_in, sampler, g).x;
  const float hxy = 2.f * h_center;

  const float exponent = 1.f + s;

  // Integration on half domain R x R-
  for (int i = -ir; i <= ir; ++i)
  {
    for (int j = -ir; j < 0; ++j)
    {
      float dist2 = (float)(i * i + j * j);
      float d = native_powr(dist2, exponent);

      int xp = clamp(g.x + i, 0, nx - 1);
      int yp = max(g.y + j, 0);

      int xn = clamp(g.x - i, 0, nx - 1);
      int yn = min(g.y - j, ny - 1);

      float hp = read_imagef(img_in, sampler, (int2)(xp, yp)).x;
      float hn = read_imagef(img_in, sampler, (int2)(xn, yn)).x;

      sum += (hp - hxy + hn) / d;
    }
  }

  // Integration on R- x {0}
  for (int i = -ir; i < 0; ++i)
  {
    float dist2 = (float)(i * i);
    float d = native_powr(dist2, exponent);

    int xp = clamp(g.x + i, 0, nx - 1);
    int xn = clamp(g.x - i, 0, nx - 1);

    float hp = read_imagef(img_in, sampler, (int2)(xp, g.y)).x;
    float hn = read_imagef(img_in, sampler, (int2)(xn, g.y)).x;

    sum += (hp - hxy + hn) / d;
  }

  write_imagef(img_out, g, sum);
}
)""
