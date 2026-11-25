R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel local_max_square(read_only image2d_t  in,
                             write_only image2d_t out,
                             const int            nx,
                             const int            ny,
                             const int            ir,
                             const int            iter)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float vmax = -FLT_MAX;

  if (iter == 0)
    for (int k = -ir; k < ir + 1; k++)
      vmax = max(vmax, read_imagef(in, sampler, (int2)(g.x + k, g.y)).x);
  else
    for (int k = -ir; k < ir + 1; k++)
      vmax = max(vmax, read_imagef(in, sampler, (int2)(g.x, g.y + k)).x);

  write_imagef(out, g, vmax);
}
)""
