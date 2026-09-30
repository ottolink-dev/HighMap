R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel topographic_position_index(read_only image2d_t  img_in,
                                       write_only image2d_t img_out,
                                       const int            nx,
                                       const int            ny,
                                       const int            ir)
{
  // Adapted from terrain-descriptors
  // Original author: oargudo
  // License: MIT (see THIRD_PARTY_LICENSES.md)

  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float s = 0.f;
  int   count = 0;

  for (int i = g.x - ir; i < g.x + ir + 1; ++i)
    for (int j = g.y - ir; j < g.y + ir + 1; ++j)
      if (is_within_radius_and_inside(i, j, g.x, g.y, ir, nx, ny))
      {
        float v = read_imagef(img_in, sampler, (int2)(i, j)).x;
        s += v;
        count++;
      }

  float val = read_imagef(img_in, sampler, g).x - s / (float)(count);

  write_imagef(img_out, g, val);
}
)""
