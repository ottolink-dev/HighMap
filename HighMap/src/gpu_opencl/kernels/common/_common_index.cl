R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#define TGET(tex, i, j) read_imagef(tex, sampler, (int2)((i), (j))).x
#define TSET(tex, i, j, v) write_imagef(tex, (int2)((i), (j)), (v))

float2 g_to_xy(const int2   g,
               const int    nx,
               const int    ny,
               const float  kx,
               const float  ky,
               const float  dx,
               const float  dy,
               const float4 bbox)
{
  float x = (float)g.x / (float)nx;
  float y = (float)g.y / (float)ny;

  x = kx * (x * (bbox.y - bbox.x) + bbox.x) + kx * dx;
  y = ky * (y * (bbox.w - bbox.z) + bbox.z) + ky * dy;

  return (float2)(x, y);
}

float2 g_to_xy_pixel_centered(const int2 g, const int nx, const int ny)
{
  float dx = 1.f / (float)nx;
  float dy = 1.f / (float)ny;

  float x = 0.5f * dx + (float)g.x * dx;
  float y = 0.5f * dy + (float)g.y * dy;
  return (float2)(x, y);
}

float2 g_to_xy_pixel_centered_bbox(const int2   g,
                                   const int    nx,
                                   const int    ny,
                                   const float4 bbox)
{
  float dx = 1.f / (float)nx;
  float dy = 1.f / (float)ny;

  float x = 0.5f * dx + (float)g.x * dx;
  float y = 0.5f * dy + (float)g.y * dy;

  x = x * (bbox.y - bbox.x) + bbox.x;
  y = y * (bbox.w - bbox.z) + bbox.z;

  return (float2)(x, y);
}

int2 xy_to_g(const float2 xy, const int nx, const int ny)
{
  return (int2)((int)(xy.x * (nx - 1.f)), (int)(xy.y * (ny - 1.f)));
}

bool is_inside(const int i, const int j, const int nx, const int ny)
{
  return i >= 0 && i < nx && j >= 0 && j < ny;
}

bool is_inside_gap(const int i,
                   const int j,
                   const int nx,
                   const int ny,
                   const int gap)
{
  return i >= gap && i < nx - gap && j >= gap && j < ny - gap;
}

bool is_within_radius_and_inside(const int i,
                                 const int j,
                                 const int ic,
                                 const int jc,
                                 const int ir,
                                 const int nx,
                                 const int ny)
{
  if (!is_inside(i, j, nx, ny)) return false;

  int id = i - ic;
  int jd = j - jc;

  return (id * id + jd * jd <= ir * ir);
}

int linear_index(const int i, const int j, const int nx)
{
  return j * nx + i;
}

int linear_index_g(const int2 g, const int nx)
{
  return g.y * nx + g.x;
}

void update_interp_param(float2 pos, int *i, int *j, float *u, float *v)
{
  *i = (int)pos.x;
  *j = (int)pos.y;
  *u = pos.x - *i;
  *v = pos.y - *j;
}

inline int apply_boundaries(global float *z, int gx, int gy, int nx, int ny)
{
  int index = linear_index(gx, gy, nx);

  if (gx == 0)
  {
    z[index] = z[linear_index(1, gy, nx)];
    return 1;
  }
  if (gx == nx - 1)
  {
    z[index] = z[linear_index(nx - 2, gy, nx)];
    return 1;
  }
  if (gy == 0)
  {
    z[index] = z[linear_index(gx, 1, nx)];
    return 1;
  }
  if (gy == ny - 1)
  {
    z[index] = z[linear_index(gx, ny - 2, nx)];
    return 1;
  }

  return 0; // not a boundary
}

inline int apply_boundaries_io(global const float *z_in,
                               global float       *z_out,
                               int                 gx,
                               int                 gy,
                               int                 nx,
                               int                 ny)
{
  int index = linear_index(gx, gy, nx);

  if (gx == 0)
  {
    z_out[index] = z_in[linear_index(1, gy, nx)];
    return 1;
  }
  if (gx == nx - 1)
  {
    z_out[index] = z_in[linear_index(nx - 2, gy, nx)];
    return 1;
  }
  if (gy == 0)
  {
    z_out[index] = z_in[linear_index(gx, 1, nx)];
    return 1;
  }
  if (gy == ny - 1)
  {
    z_out[index] = z_in[linear_index(gx, ny - 2, nx)];
    return 1;
  }

  return 0; // not a boundary
}

inline int apply_boundaries_buffer(global float *z,
                                   int           gx,
                                   int           gy,
                                   int           nx,
                                   int           ny,
                                   int           b)
{
  // inside valid interior → nothing to do
  if (gx >= b && gx < nx - b && gy >= b && gy < ny - b) return 0;

  // clamp to nearest interior cell
  int cx = clamp(gx, b, nx - b - 1);
  int cy = clamp(gy, b, ny - b - 1);

  int idx = linear_index(gx, gy, nx);
  int idxc = linear_index(cx, cy, nx);

  z[idx] = z[idxc];
  return 1;
}

inline int apply_boundaries_buffer_io(global const float *z_in,
                                      global float       *z_out,
                                      int                 gx,
                                      int                 gy,
                                      int                 nx,
                                      int                 ny,
                                      int                 b)
{
  // inside valid interior → nothing to do
  if (gx >= b && gx < nx - b && gy >= b && gy < ny - b) return 0;

  // clamp to nearest interior cell
  int cx = clamp(gx, b, nx - b - 1);
  int cy = clamp(gy, b, ny - b - 1);

  int idx = linear_index(gx, gy, nx);
  int idxc = linear_index(cx, cy, nx);

  z_out[idx] = z_in[idxc];
  return 1;
}
)""
