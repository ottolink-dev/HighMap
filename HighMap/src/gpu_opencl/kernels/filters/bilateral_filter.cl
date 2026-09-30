R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel bilateral_filter(read_only image2d_t  array,
                             read_only image2d_t  kernel2d,
                             read_only image2d_t  kernel1d,
                             write_only image2d_t out,
                             const int            nx,
                             const int            ny,
                             const int            k2dx,
                             const int            k2dy,
                             const int            k1dx,
                             const float          kernel1d_value_scaling)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const sampler_t sampler_itp = CLK_NORMALIZED_COORDS_TRUE |
                                CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  const int i = g.x;
  const int j = g.y;

  const int   dx = (int)(0.5f * k2dx);
  const int   dy = (int)(0.5f * k2dy);
  const float val_center = TGET(array, i, j);

  float new_val = 0.f;
  float weight = 0.f;

  // brute-force filter
  for (int r = 0; r < k2dx; r++)
    for (int s = 0; s < k2dy; s++)
    {
      float val_current = TGET(array, i + r - dx, j + s - dy);
      float diff = val_current - val_center;

      // stretch assuming diff in [-1..1]
      diff = clamp(diff * 0.5f * k1dx * kernel1d_value_scaling,
                   -(float)k1dx,
                   (float)k1dx);

      float weight2d = TGET(kernel2d, r, s);
      float weight1d = read_imagef(kernel1d, sampler_itp, (float2)(diff, 0.f))
                           .x;

      new_val += val_current * weight2d * weight1d;
      weight += weight2d * weight1d;
    }

  if (weight > 1e-6f)
    new_val /= weight;
  else
    new_val = 0.f;

  TSET(out, i, j, new_val);
}
)""
