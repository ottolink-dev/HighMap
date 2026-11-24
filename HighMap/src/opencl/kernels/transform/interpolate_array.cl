R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#define MAX_LAGRANGE_SIZE 16

void kernel interpolate_array_bicubic(read_only image2d_t  source,
                                      write_only image2d_t target,
                                      int                  nx_s,
                                      int                  ny_s,
                                      int                  nx_t,
                                      int                  ny_t,
                                      float4               bbox_t)
{
  // kernel distributed over the target pixels (ie. nx_t * ny_t calculations)
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx_t || g.y >= ny_t) return;

  // --- find nearest point in source from interpolation position

  float dx_s = 1.f / (float)nx_s;
  float dy_s = 1.f / (float)ny_s;

  float dx_t = 1.f / (float)nx_t;
  float dy_t = 1.f / (float)ny_t;

  // pixel-centered coordinates
  float x = (0.5f + (float)g.x) * dx_t;
  float y = (0.5f + (float)g.y) * dy_t;

  // scale and shift according to the bounding box
  x = x * (bbox_t.y - bbox_t.x) + bbox_t.x;
  y = y * (bbox_t.w - bbox_t.z) + bbox_t.z;

  // corresponding index for source array
  x = (x - 0.5f * dx_s) / dx_s;
  y = (y - 0.5f * dy_s) / dy_s;

  int2 g_s = (int2)((int)x, (int)y);

  float u = x - g_s.x;
  float v = y - g_s.y;

  // --- interpolate

  float arr[4][4];

  // get the 4x4 surrounding grid points
  for (int n = -1; n <= 2; ++n)
    for (int m = -1; m <= 2; ++m)
    {
      int ip = clamp(g_s.x + m, 0, nx_s - 1);
      int jp = clamp(g_s.y + n, 0, ny_s - 1);

      arr[m + 1][n + 1] = read_imagef(source, (int2)(ip, jp)).x;
    }

  // interpolate in the x direction
  float col_results[4];
  for (int k = 0; k < 4; ++k)
    col_results[k] = cubic_interp_vec(arr[k], v);

  // interpolate in the y direction
  float val = cubic_interp_vec(col_results, u);

  write_imagef(target, g, val);
}

void kernel interpolate_array_bilinear(read_only image2d_t  source,
                                       write_only image2d_t target,
                                       int                  nx_t,
                                       int                  ny_t,
                                       float4               bbox_t)
{
  // kernel distributed over the target pixels (ie. nx_t * ny_t calculations)
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx_t || g.y >= ny_t) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_TRUE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  // bounding box of the source array is assumed to be a unit square
  float2 pos_target = g_to_xy_pixel_centered_bbox(g, nx_t, ny_t, bbox_t);

  float val = read_imagef(source, sampler, pos_target).x;

  write_imagef(target, g, val);
}

void kernel interpolate_array_nearest(read_only image2d_t  source,
                                      write_only image2d_t target,
                                      int                  nx_t,
                                      int                  ny_t,
                                      float4               bbox_t)
{
  // kernel distributed over the target pixels (ie. nx_t * ny_t calculations)
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx_t || g.y >= ny_t) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_TRUE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  // bounding box of the source array is assumed to be a unit square
  float2 pos_target = g_to_xy_pixel_centered_bbox(g, nx_t, ny_t, bbox_t);

  float val = read_imagef(source, sampler, pos_target).x;

  write_imagef(target, g, val);
}

//
void kernel interpolate_array_lagrange(read_only image2d_t  source,
                                       write_only image2d_t target,
                                       int                  nx_s,
                                       int                  ny_s,
                                       int                  nx_t,
                                       int                  ny_t,
                                       int                  order)
{
  // kernel distributed over the target pixels (ie. nx_t * ny_t calculations)
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx_t || g.y >= ny_t) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_TRUE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  float2 pos_target = g_to_xy_pixel_centered(g, nx_t, ny_t);

  // --- find nearest point in source from interpolation position

  float dx_s = 1.f / (float)nx_s;
  float dy_s = 1.f / (float)ny_s;

  float dx_t = 1.f / (float)nx_t;
  float dy_t = 1.f / (float)ny_t;

  // pixel-centered coordinates
  float x = (0.5f + (float)g.x) * dx_t;
  float y = (0.5f + (float)g.y) * dy_t;

  // node-centered coordinates with index scaling (0, 1, 2...)
  x = (x - 0.5f * dx_s) / dx_s;
  y = (y - 0.5f * dy_s) / dy_s;
  float2 pos = (float2)(x, y);

  // corresponding index for source array
  int2 g_s = (int2)((int)pos.x, (int)pos.y);

  // --- compute Lagrangian coefficients

  float lx[MAX_LAGRANGE_SIZE];
  float ly[MAX_LAGRANGE_SIZE];
  int   nk = max(1, (int)(0.5f * (order - 1)));

  // x-direction
  for (int i = 0; i < order; ++i)
  {
    lx[i] = 1.f;
    float xi = (g_s.x - nk + i);

    for (int k = 0; k < order; ++k)
      if (k != i)
      {
        float xk = (g_s.x - nk + k);
        lx[i] *= (pos.x - xk) / (xi - xk);
      }
  }

  // y-direction
  for (int j = 0; j < order; ++j)
  {
    ly[j] = 1.f;
    float yj = (g_s.y - nk + j);

    for (int k = 0; k < order; ++k)
      if (k != j)
      {
        float yk = (g_s.y - nk + k);
        ly[j] *= (pos.y - yk) / (yj - yk);
      }
  }

  // --- interpolate

  float val = 0.f;

  for (int i = 0; i < order; ++i)
    for (int j = 0; j < order; ++j)
    {
      int ip = clamp(g_s.x - nk + i, 0, nx_s - 1);
      int jp = clamp(g_s.y - nk + j, 0, ny_s - 1);

      val += lx[i] * ly[j] * read_imagef(source, (int2)(ip, jp)).x;
    }

  write_imagef(target, g, val);
}
)""
