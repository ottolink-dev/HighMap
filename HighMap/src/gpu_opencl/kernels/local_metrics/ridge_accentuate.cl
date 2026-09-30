R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel ridge_accentuate(read_only image2d_t  array,
                             read_only image2d_t  array_f,
                             write_only image2d_t out,
                             const int            nx,
                             const int            ny,
                             const float          strength)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler_nearest = CLK_NORMALIZED_COORDS_FALSE |
                                    CLK_ADDRESS_CLAMP_TO_EDGE |
                                    CLK_FILTER_NEAREST;

  const sampler_t sampler_interp = CLK_NORMALIZED_COORDS_FALSE |
                                   CLK_ADDRESS_CLAMP_TO_EDGE |
                                   CLK_FILTER_LINEAR;

  float scale_x = (float)(nx - 1);
  float scale_y = (float)(ny - 1);

  // normal vector from central differences of smoothed terrain scaled to domain
  // units
  float dx = 0.5f *
             (read_imagef(array_f, sampler_nearest, (int2)(g.x + 1, g.y)).x -
              read_imagef(array_f, sampler_nearest, (int2)(g.x - 1, g.y)).x) *
             scale_x;
  float dy = 0.5f *
             (read_imagef(array_f, sampler_nearest, (int2)(g.x, g.y + 1)).x -
              read_imagef(array_f, sampler_nearest, (int2)(g.x, g.y - 1)).x) *
             scale_y;

  float3 normal = {-dx, -dy, 1.f};
  normal /= length(normal);

  // displacement along horizontal normal direction
  // in OpenCL linear filtering with unnormalized coordinates, pixel centers are
  // at g + 0.5
  float2 pos = {(float)g.x + 0.5f + strength * scale_x * normal.x,
                (float)g.y + 0.5f + strength * scale_y * normal.y};
  float  val = read_imagef(array, sampler_interp, pos).x;

  write_imagef(out, g, val);
}
)""
