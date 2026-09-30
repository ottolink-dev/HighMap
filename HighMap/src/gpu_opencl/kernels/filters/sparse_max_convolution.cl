R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel sparse_max_convolution(global float       *array,
                                   read_only image2d_t weights,
                                   global float       *out,
                                   const int           nx,
                                   const int           ny,
                                   const int           nkx,
                                   const int           nky,
                                   const int           scale_kernel_with_value,
                                   const float         k_smoothmax)
{
  int2 g = (int2)(get_global_id(0), get_global_id(1));
  if (g.x >= nx || g.y >= ny) return;

  int   index = linear_index(g.x, g.y, nx);
  float val = array[index];

  if (val == 0.f) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  // local kernel size
  int max_irx = nkx / 2;
  int max_iry = nky / 2;

  // scale kernel radius according to val
  int irx = nkx / 2;
  int iry = nky / 2;

  if (scale_kernel_with_value)
  {
    irx = max(1, (int)round((nkx / 2.0f) * val));
    iry = max(1, (int)round((nky / 2.0f) * val));
  }

  for (int dy = -iry; dy <= iry; ++dy)
    for (int dx = -irx; dx <= irx; ++dx)
    {
      int i = g.x + dx;
      int j = g.y + dy;

      if (i < 0 || i >= nx || j < 0 || j >= ny) continue;

      // map the local kernel coordinates back into the original
      // weight texture coordinates.
      float u = (dx + irx) / (float)(2 * irx);
      float v = (dy + iry) / (float)(2 * iry);

      float xw = u * (nkx - 1) + 0.5f;
      float yw = v * (nky - 1) + 0.5f;

      float w = read_imagef(weights, sampler, (float2)(xw, yw)).x * val;

      atomic_smoothmax_float(&out[linear_index(i, j, nx)], w, k_smoothmax);
    }
}

void kernel sparse_max_convolution_with_size(global float       *array_values,
                                             global float       *array_sizes,
                                             read_only image2d_t weights,
                                             global float       *out,
                                             const int           nx,
                                             const int           ny,
                                             const int           nkx,
                                             const int           nky,
                                             const float         k_smoothmax)
{
  int2 g = (int2)(get_global_id(0), get_global_id(1));
  if (g.x >= nx || g.y >= ny) return;

  int   index = linear_index(g.x, g.y, nx);
  float val = array_values[index];
  float size = array_sizes[index];

  if (val == 0.f) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  // local kernel size
  int max_irx = nkx / 2;
  int max_iry = nky / 2;

  // scale kernel radius according to val
  int irx = nkx / 2;
  int iry = nky / 2;

  irx = max(1, (int)round((nkx / 2.0f) * size));
  iry = max(1, (int)round((nky / 2.0f) * size));

  for (int dy = -iry; dy <= iry; ++dy)
    for (int dx = -irx; dx <= irx; ++dx)
    {
      int i = g.x + dx;
      int j = g.y + dy;

      if (i < 0 || i >= nx || j < 0 || j >= ny) continue;

      // map the local kernel coordinates back into the original
      // weight texture coordinates.
      float u = (dx + irx) / (float)(2 * irx);
      float v = (dy + iry) / (float)(2 * iry);

      float xw = u * (nkx - 1) + 0.5f;
      float yw = v * (nky - 1) + 0.5f;

      float w = read_imagef(weights, sampler, (float2)(xw, yw)).x * val;

      atomic_smoothmax_float(&out[linear_index(i, j, nx)], w, k_smoothmax);
    }
}
)""
