R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel smooth_cpulse(read_only image2d_t  in,
                          read_only image2d_t  w,
                          write_only image2d_t out,
                          const int            nx,
                          const int            ny,
                          const int            ir,
                          const int            pass_nb)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float val = 0.f;

  if (pass_nb == 0)
  {
    for (int k = -ir; k < ir + 1; k++)
    {
      float f = read_imagef(in, sampler, (int2)(g.x + k, g.y)).x;
      f *= read_imagef(w, sampler, (int2)(k + ir, 0)).x;
      val += f;
    }
  }
  else
  {
    for (int k = -ir; k < ir + 1; k++)
    {
      float f = read_imagef(in, sampler, (int2)(g.x, g.y + k)).x;
      f *= read_imagef(w, sampler, (int2)(k + ir, 0)).x;
      val += f;
    }
  }

  write_imagef(out, (int2)(g.x, g.y), val);
}
)""
