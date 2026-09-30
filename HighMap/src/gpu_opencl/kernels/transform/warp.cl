R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel warp_xy(read_only image2d_t  in,
                    read_only image2d_t  dx,
                    read_only image2d_t  dy,
                    write_only image2d_t out,
                    int                  nx,
                    int                  ny)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler_in = CLK_NORMALIZED_COORDS_TRUE |
                               CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

  const sampler_t sampler_dxy = CLK_NORMALIZED_COORDS_FALSE |
                                CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float vdx = read_imagef(dx, sampler_dxy, g).x;
  float vdy = read_imagef(dy, sampler_dxy, g).x;

  float2 pos = {(float)g.x / (nx - 1.f) + vdx, (float)g.y / (ny - 1.f) + vdy};

  float val = read_imagef(in, sampler_in, pos).x;

  write_imagef(out, g, val);
}

void kernel warp_x(read_only image2d_t  in,
                   read_only image2d_t  dx,
                   write_only image2d_t out,
                   int                  nx,
                   int                  ny)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler_in = CLK_NORMALIZED_COORDS_TRUE |
                               CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

  const sampler_t sampler_dxy = CLK_NORMALIZED_COORDS_FALSE |
                                CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float vdx = read_imagef(dx, sampler_dxy, g).x;

  float2 pos = {(float)g.x / (nx - 1.f) + vdx, (float)g.y / (ny - 1.f)};

  float val = read_imagef(in, sampler_in, pos).x;

  write_imagef(out, g, val);
}

void kernel warp_y(read_only image2d_t  in,
                   read_only image2d_t  dy,
                   write_only image2d_t out,
                   int                  nx,
                   int                  ny)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler_in = CLK_NORMALIZED_COORDS_TRUE |
                               CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

  const sampler_t sampler_dxy = CLK_NORMALIZED_COORDS_FALSE |
                                CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float vdy = read_imagef(dy, sampler_dxy, g).x;

  float2 pos = {(float)g.x / (nx - 1.f), (float)g.y / (ny - 1.f) + vdy};

  float val = read_imagef(in, sampler_in, pos).x;

  write_imagef(out, g, val);
}
)""
