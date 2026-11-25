R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel normal_displacement(read_only image2d_t  array,
                                read_only image2d_t  array_f,
                                write_only image2d_t out,
                                const int            nx,
                                const int            ny,
                                const float          amount)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler_nearest = CLK_NORMALIZED_COORDS_FALSE |
                                    CLK_ADDRESS_CLAMP_TO_EDGE |
                                    CLK_FILTER_NEAREST;

  const sampler_t sampler_interp = CLK_NORMALIZED_COORDS_TRUE |
                                   CLK_ADDRESS_MIRRORED_REPEAT |
                                   CLK_FILTER_LINEAR;

  // shape factor
  float cx = 0.5f * (nx - 1.f);
  float cy = 0.5f * (ny - 1.f);
  float x = (g.x - cx) / cx;
  float y = (g.y - cy) / cy;
  float r2 = x * x + y * y;
  float factor = r2 < 1.f ? (1.f - r2) * (1.f - r2) : 0.f;

  // normal vector
  float dx = read_imagef(array_f, sampler_nearest, (int2)(g.x + 1, g.y)).x -
             read_imagef(array_f, sampler_nearest, (int2)(g.x - 1, g.y)).x;
  float dy = read_imagef(array_f, sampler_nearest, (int2)(g.x, g.y + 1)).x -
             read_imagef(array_f, sampler_nearest, (int2)(g.x, g.y - 1)).x;

  float3 normal = {-0.5f * dx, -0.5f * dy, 1.f};
  normal /= length(normal);
  normal *= -amount * factor;

  // warp
  float  val0 = read_imagef(array, sampler_nearest, g).x;
  float2 pos = {(float)g.x / (nx - 1.f) + normal.x,
                (float)g.y / (ny - 1.f) + normal.y};
  float  val = read_imagef(array, sampler_interp, pos).x;

  // write output
  write_imagef(out, g, val);
}
)""
