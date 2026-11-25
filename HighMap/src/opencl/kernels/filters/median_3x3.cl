R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel median_3x3(read_only image2d_t  img_in,
                       write_only image2d_t img_out,
                       int                  nx,
                       int                  ny)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float neighborhood[9];

  neighborhood[0] = read_imagef(img_in, sampler, (int2)(g.x - 1, g.y - 1)).x;
  neighborhood[1] = read_imagef(img_in, sampler, (int2)(g.x, g.y - 1)).x;
  neighborhood[2] = read_imagef(img_in, sampler, (int2)(g.x + 1, g.y - 1)).x;
  neighborhood[3] = read_imagef(img_in, sampler, (int2)(g.x - 1, g.y)).x;
  neighborhood[4] = read_imagef(img_in, sampler, (int2)(g.x, g.y)).x;
  neighborhood[5] = read_imagef(img_in, sampler, (int2)(g.x + 1, g.y)).x;
  neighborhood[6] = read_imagef(img_in, sampler, (int2)(g.x - 1, g.y + 1)).x;
  neighborhood[7] = read_imagef(img_in, sampler, (int2)(g.x, g.y + 1)).x;
  neighborhood[8] = read_imagef(img_in, sampler, (int2)(g.x + 1, g.y + 1)).x;
  insertion_sort(neighborhood, 9);

  const float val = neighborhood[4];

  write_imagef(img_out, (int2)(g.x, g.y), val);
}
)""
