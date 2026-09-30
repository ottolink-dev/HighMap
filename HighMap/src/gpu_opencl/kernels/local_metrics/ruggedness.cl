R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel ruggedness(read_only image2d_t  array,
                       write_only image2d_t out,
                       const int            nx,
                       const int            ny,
                       const int            ir)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float val = 0.f;
  float center = read_imagef(array, sampler, g).x;

  for (int i = g.x - ir; i < g.x + ir + 1; i++)
    for (int j = g.y - ir; j < g.y + ir + 1; j++)
      if (is_within_radius_and_inside(i, j, g.x, g.y, ir, nx, ny))
      {
        float delta = center - read_imagef(array, sampler, (int2)(i, j)).x;
        val += delta * delta;
      }

  write_imagef(out, g, sqrt(val));
}
)""
