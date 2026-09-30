R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel expand(read_only image2d_t  z,
                   read_only image2d_t  w,
                   write_only image2d_t out,
                   const int            nx,
                   const int            ny,
                   const int            kx,
                   const int            ky)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float val = 0.f;
  int   dx = (int)(0.5f * kx);
  int   dy = (int)(0.5f * ky);

  for (int k = 0; k < kx; k++)
    for (int r = 0; r < ky; r++)
    {
      float f = read_imagef(z, sampler, (int2)(g.x + k - dx, g.y + r - dy)).x;
      f *= read_imagef(w, sampler, (int2)(k, r)).x;
      val = max(val, f);
    }

  write_imagef(out, (int2)(g.x, g.y), val);
}
)""
