R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#define ROCK_STATUS_ACTIVE 0
#define ROCK_STATUS_RESTING 1
#define ROCK_STATUS_OUT_OF_BOUNDS 2

// --- Bilinear Sampling Helper

inline float rock_helper_sample_height(global const float *z,
                                       int                 i,
                                       int                 j,
                                       int                 nx,
                                       int                 ny,
                                       float               u,
                                       float               v)
{
  if (i >= nx - 1) i = nx - 2;
  if (j >= ny - 1) j = ny - 2;
  if (i < 0) i = 0;
  if (j < 0) j = 0;

  float f00 = z[linear_index(i, j, nx)];
  float f10 = z[linear_index(i + 1, j, nx)];
  float f01 = z[linear_index(i, j + 1, nx)];
  float f11 = z[linear_index(i + 1, j + 1, nx)];

  float a10 = f10 - f00;
  float a01 = f01 - f00;
  float a11 = f11 - f10 - f01 + f00;

  return f00 + a10 * u + a01 * v + a11 * u * v;
}

inline float2 rock_helper_sample_gradient(global const float *z,
                                          int                 i,
                                          int                 j,
                                          int                 nx,
                                          int                 ny,
                                          float               u,
                                          float               v,
                                          float4              bbox)
{
  if (i >= nx - 1) i = nx - 2;
  if (j >= ny - 1) j = ny - 2;
  if (i < 0) i = 0;
  if (j < 0) j = 0;

  float f00 = z[linear_index(i, j, nx)];
  float f10 = z[linear_index(i + 1, j, nx)];
  float f01 = z[linear_index(i, j + 1, nx)];
  float f11 = z[linear_index(i + 1, j + 1, nx)];

  float a11 = f11 - f10 - f01 + f00;
  float df_du = (f10 - f00) + a11 * v;
  float df_dv = (f01 - f00) + a11 * u;

  // Scale slope by domain bounding box (slope = 1 => 1 unit elevation rise
  // across domain)
  float dom_w = bbox.y - bbox.x;
  float dom_h = bbox.w - bbox.z;
  float d_gx = (nx > 1 && dom_w > 1e-6f) ? ((float)(nx - 1) / dom_w) : 1.0f;
  float d_gy = (ny > 1 && dom_h > 1e-6f) ? ((float)(ny - 1) / dom_h) : 1.0f;

  return (float2)(df_du * d_gx, df_dv * d_gy);
}

// ============================================================================
//  Rock Rolling Simulation Kernel
// ============================================================================

kernel void rock_simulate_physics(global float4       *pos_rad,
                                  global float4       *vel_mass,
                                  global uint2        *status_class,
                                  global const float4 *init_pos_rad,
                                  global const float  *z,
                                  global const float  *friction_map,
                                  int                  num_rocks,
                                  int                  nx,
                                  int                  ny,
                                  float4               bbox,
                                  float                dt,
                                  float                gravity,
                                  float                soil_friction,
                                  float                rolling_resistance,
                                  float                inter_rock_restitution,
                                  float                min_velocity,
                                  int                  has_friction_map,
                                  int                  respawn_out_of_bounds,
                                  uint                 seed,
                                  int                  sub_steps)
{
  int rock_idx = get_global_id(0);
  if (rock_idx >= num_rocks) return;

  uint2 sc = status_class[rock_idx];
  if (sc.x == ROCK_STATUS_OUT_OF_BOUNDS) return;

  float4 p = pos_rad[rock_idx];
  float4 vm = vel_mass[rock_idx];

  float2 pos = (float2)(p.x, p.y);
  float  radius = p.w;
  float2 vel = (float2)(vm.x, vm.y);

  float dom_w = bbox.y - bbox.x;
  float dom_h = bbox.w - bbox.z;

  uint rng_state = wang_hash((uint)rock_idx + seed * 1999u + 1u);

  for (int step = 0; step < sub_steps; ++step)
  {
    // Check if rock is out of domain bounds
    if (pos.x < bbox.x || pos.x > bbox.y || pos.y < bbox.z || pos.y > bbox.w)
    {
      if (respawn_out_of_bounds != 0)
      {
        rng_state = wang_hash(rng_state + (uint)step * 104729u);
        uint   random_rock_idx = rng_state % (uint)num_rocks;
        float4 ip = init_pos_rad[random_rock_idx];
        pos = (float2)(ip.x, ip.y);
        vel = (float2)(0.0f, 0.0f);
        sc.x = ROCK_STATUS_ACTIVE;
      }
      else
      {
        sc.x = ROCK_STATUS_OUT_OF_BOUNDS;
        break;
      }
    }

    // Sample terrain gradient and height
    float u_norm = (pos.x - bbox.x) / dom_w;
    float v_norm = (pos.y - bbox.z) / dom_h;
    float grid_x = clamp(u_norm * (float)(nx - 1), 0.0f, (float)(nx - 1));
    float grid_y = clamp(v_norm * (float)(ny - 1), 0.0f, (float)(ny - 1));

    int   gi, gj;
    float gu, gv;
    update_interp_param((float2)(grid_x, grid_y), &gi, &gj, &gu, &gv);

    float2 grad_z =
        rock_helper_sample_gradient(z, gi, gj, nx, ny, gu, gv, bbox);

    // Downhill rolling acceleration along slope surface
    float slope_sq = dot(grad_z, grad_z);
    float inv_denom = 1.0f / (1.0f + slope_sq);
    float cos_theta = sqrt(inv_denom);
    float g_n = gravity * cos_theta;

    // Rolling sphere linear acceleration factor = 5/7
    float2 a_grav = -(5.0f / 7.0f) * gravity * grad_z * inv_denom;
    float  a_grav_mag = length(a_grav);

    // Effective friction coefficient (soil friction + rolling resistance)
    float mu_soil = soil_friction;
    if (has_friction_map != 0)
    {
      mu_soil = rock_helper_sample_height(friction_map, gi, gj, nx, ny, gu, gv);
    }
    float mu_eff = mu_soil + rolling_resistance;

    float  v_speed = length(vel);
    float2 accel = (float2)(0.0f, 0.0f);

    if (v_speed > 1e-4f)
    {
      float2 v_dir = vel / v_speed;
      float  f_mag = mu_eff * g_n;
      float2 f_friction = -min(f_mag, v_speed / dt) * v_dir;
      accel = a_grav + f_friction;
    }
    else
    {
      // Static rolling resistance threshold
      if (a_grav_mag <= mu_soil * g_n)
      {
        vel = (float2)(0.0f, 0.0f);
        accel = (float2)(0.0f, 0.0f);
      }
      else
      {
        accel = a_grav - (mu_soil * g_n) * (a_grav / a_grav_mag);
      }
    }

    // Explicit Euler integration for rolling motion
    vel += accel * dt;
    pos += vel * dt;

    // --- Inter-Rock Collisions & Overlap Resolution
    // Resolve pairwise contact / overlap with other active/resting rocks
    for (int other_idx = 0; other_idx < num_rocks; ++other_idx)
    {
      if (other_idx == rock_idx) continue;

      uint2 other_sc = status_class[other_idx];
      if (other_sc.x == ROCK_STATUS_OUT_OF_BOUNDS) continue;

      float4 other_p = pos_rad[other_idx];
      float2 other_pos = (float2)(other_p.x, other_p.y);
      float  other_radius = other_p.w;

      float2 delta = pos - other_pos;
      float  dist_sq = dot(delta, delta);
      float  min_dist = radius + other_radius;
      float  min_dist_sq = min_dist * min_dist;

      if (dist_sq < min_dist_sq && dist_sq > 1e-12f)
      {
        float  dist = sqrt(dist_sq);
        float2 n = delta / dist;
        float  overlap = min_dist - dist;

        // Position relaxation: push apart so rocks do not overlap
        pos += 0.5f * overlap * n;

        // Relative velocity collision response (impulse with restitution)
        float4 other_vm = vel_mass[other_idx];
        float2 other_vel = (float2)(other_vm.x, other_vm.y);
        float2 rel_vel = vel - other_vel;
        float  v_norm = dot(rel_vel, n);

        if (v_norm < 0.0f)
        {
          float m1 = vm.w;
          float m2 = other_vm.w;
          float inv_m1 = (m1 > 1e-6f) ? (1.0f / m1) : 1.0f;
          float inv_m2 = (m2 > 1e-6f) ? (1.0f / m2) : 1.0f;

          float impulse = -(1.0f + inter_rock_restitution) * v_norm /
                          (inv_m1 + inv_m2);
          vel += (impulse * inv_m1) * n;
        }
      }
    }

    // Check resting state on gentle slope / flat
    v_speed = length(vel);
    if (v_speed < min_velocity && a_grav_mag <= mu_soil * g_n)
    {
      vel = (float2)(0.0f, 0.0f);
      sc.x = ROCK_STATUS_RESTING;
    }
    else
    {
      sc.x = ROCK_STATUS_ACTIVE;
    }
  }

  // Sample final terrain elevation for continuous surface contact
  float u_norm = (pos.x - bbox.x) / dom_w;
  float v_norm = (pos.y - bbox.z) / dom_h;
  float grid_x = clamp(u_norm * (float)(nx - 1), 0.0f, (float)(nx - 1));
  float grid_y = clamp(v_norm * (float)(ny - 1), 0.0f, (float)(ny - 1));
  int   gi, gj;
  float gu, gv;
  update_interp_param((float2)(grid_x, grid_y), &gi, &gj, &gu, &gv);
  float z_final = rock_helper_sample_height(z, gi, gj, nx, ny, gu, gv) + radius;

  pos_rad[rock_idx] = (float4)(pos.x, pos.y, z_final, radius);
  vel_mass[rock_idx] = (float4)(vel.x, vel.y, 0.0f, vm.w);
  status_class[rock_idx] = sc;
}
)""
