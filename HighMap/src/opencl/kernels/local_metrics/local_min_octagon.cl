R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel local_min_octagon(read_only image2d_t  in,
                              write_only image2d_t out,
                              const int            nx,
                              const int            ny,
                              const int            radius,
                              const int            iter)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float vmin = FLT_MAX;

  if (iter == 0) // horizontal pass
    for (int k = -radius; k <= radius; k++)
      vmin = min(vmin, read_imagef(in, sampler, (int2)(g.x + k, g.y)).x);
  else if (iter == 1) // vertical pass
    for (int k = -radius; k <= radius; k++)
      vmin = min(vmin, read_imagef(in, sampler, (int2)(g.x, g.y + k)).x);
  else if (iter == 2) // main diagonal (+1, +1) pass
    for (int k = -radius; k <= radius; k++)
      vmin = min(vmin, read_imagef(in, sampler, (int2)(g.x + k, g.y + k)).x);
  else if (iter == 3) // anti diagonal (+1, -1) pass
    for (int k = -radius; k <= radius; k++)
      vmin = min(vmin, read_imagef(in, sampler, (int2)(g.x + k, g.y - k)).x);

  write_imagef(out, g, vmin);
}
)""
