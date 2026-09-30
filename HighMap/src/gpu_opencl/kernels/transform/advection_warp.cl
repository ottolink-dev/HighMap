R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#define MAX_STEPS 1024 // must be >= advection_length / min(drx,dry)

void kernel advection_warp(read_only image2d_t  z,
                           read_only image2d_t  field,
                           read_only image2d_t  dx,
                           read_only image2d_t  dy,
                           read_only image2d_t  mask,
                           write_only image2d_t out,
                           const int            nx,
                           const int            ny,
                           const float          advection_length,
                           const float          value_persistence)
{
  const int2 g = {get_global_id(0), get_global_id(1)};
  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_TRUE |
                            CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

  const float drx = 1.0f / (float)nx;
  const float dry = 1.0f / (float)ny;
  const float eps = 0.f; // 1e-6f;

  const float2 pos0 = g_to_xy_pixel_centered(g, nx, ny);
  float2       pos = pos0;

  float     z_prev = read_imagef(z, sampler, pos).x;
  const int max_steps = (int)ceil(advection_length / fmin(drx, dry));
  const int limit = min(max_steps, MAX_STEPS);

  // local storage of positions
  float2 path[MAX_STEPS];
  int    path_len = 0;

  // ==========================================================
  // STEP 1 : upward pass — record all visited positions
  // ==========================================================

  float2 dir;

  for (int i = 0; i < limit; ++i)
  {
    path[path_len++] = pos; // store current position

    float2 vel = (float2)(read_imagef(dx, sampler, pos).x,
                          read_imagef(dy, sampler, pos).x);

    float len = hypot(vel.x, vel.y);
    if (len > eps)
    {
      // (NB - if too flat, keep the direction of the previous step)
      dir = vel / len;
    }

    pos.x += dir.x * drx;
    pos.y += dir.y * dry;

    if (pos.x < 0.f || pos.x > 1.0f || pos.y < 0.f || pos.y > 1.0f) break;

    float z_current = read_imagef(z, sampler, pos).x;
    if (z_current - z_prev < -1e-5f) break;

    z_prev = z_current;
  }

  // ==========================================================
  // STEP 2 : downward pass — walk backward through stored path
  // ==========================================================

  float val = 0.0f;
  if (path_len > 0) val = read_imagef(field, sampler, path[path_len - 1]).x;

  for (int i = path_len - 1; i >= 0; --i)
  {
    pos = path[i];

    float new_val = read_imagef(field, sampler, pos).x;
    float mask_v = read_imagef(mask, sampler, pos).x;
    float ratio = mask_v * value_persistence;

    val = ratio * val + (1.0f - ratio) * new_val;
  }

  write_imagef(out, g, val);
}
)""
