R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// --- HELPERS

inline int helper_wrap_coord(int c, int n)
{
  int r = c % n;
  return r < 0 ? r + n : r;
}

inline float helper_sample_wrapped_height(global const float *z, int x, int y, int nx, int ny)
{
  int wx = helper_wrap_coord(x, nx);
  int wy = helper_wrap_coord(y, ny);
  return z[linear_index(wx, wy, nx)];
}

inline bool helper_is_in_shadow(global const float *z,
                                int                 lx,
                                int                 ly,
                                float               z_land,
                                float2              wdir,
                                float               shadow_talus,
                                int                 nx,
                                int                 ny,
                                int                 periodic)
{
  // march upwind along -wdir
  float2 upwind_dir = -wdir;
  int max_steps = min(32, max(nx, ny) / 4);

  for (int step = 1; step <= max_steps; ++step)
  {
    float ux = (float)lx + upwind_dir.x * (float)step;
    float uy = (float)ly + upwind_dir.y * (float)step;

    int ix = (int)floor(ux + 0.5f);
    int iy = (int)floor(uy + 0.5f);

    float z_up = 0.f;
    if (periodic != 0)
    {
      z_up = helper_sample_wrapped_height(z, ix, iy, nx, ny);
    }
    else
    {
      if (!is_inside(ix, iy, nx, ny)) break;
      z_up = z[linear_index(ix, iy, nx)];
    }

    // if upwind terrain is high enough to cast shadow with shadow talus threshold
    float dist = (float)step;
    if (z_up - z_land > dist * shadow_talus)
    {
      return true;
    }
  }

  return false;
}

// --- KERNELS

void kernel sand_dune_step(global float       *z,
                           global const float *bedrock,
                           const int           nx,
                           const int           ny,
                           const int           nparticles,
                           const uint          seed,
                           const float         wind_angle,
                           const float         wind_speed,
                           const float         crest_speedup,
                           const float         hop_length,
                           const float         elevation_hop_factor,
                           const float         jitter_angle,
                           const float         p_sand,
                           const float         p_bare,
                           const float         shadow_talus,
                           const float         slab_height,
                           const int           max_hops,
                           const int           periodic,
                           const float         sand_inflow_rate,
                           const int           has_bedrock,
                           const int           pass_id)
{
  int gid = get_global_id(0);
  if (gid >= nparticles) return;

  uint rng = wang_hash(seed ^ (gid * 1973u + pass_id * 9277u + 11u));

  // determine launch position
  int i0 = (int)(rand(&rng) * (float)nx);
  int j0 = (int)(rand(&rng) * (float)ny);
  i0 = clamp(i0, 0, nx - 1);
  j0 = clamp(j0, 0, ny - 1);

  int idx0 = linear_index(i0, j0, nx);
  float z0 = z[idx0];

  // check bedrock constraint at launch cell
  float slab = slab_height;
  if (has_bedrock != 0)
  {
    float sand_avail = z0 - bedrock[idx0];
    if (sand_avail <= 1e-5f) return;
    slab = min(slab, sand_avail);
  }

  // atomically erode slab at launch cell
  atomic_add_float(&z[idx0], -slab);

  // calculate local crest speedup
  float z_north = (periodic != 0) ? helper_sample_wrapped_height(z, i0, j0 + 1, nx, ny)
                                  : z[linear_index(i0, min(j0 + 1, ny - 1), nx)];
  float z_south = (periodic != 0) ? helper_sample_wrapped_height(z, i0, j0 - 1, nx, ny)
                                  : z[linear_index(i0, max(j0 - 1, 0), nx)];
  float z_east  = (periodic != 0) ? helper_sample_wrapped_height(z, i0 + 1, j0, nx, ny)
                                  : z[linear_index(min(i0 + 1, nx - 1), j0, nx)];
  float z_west  = (periodic != 0) ? helper_sample_wrapped_height(z, i0 - 1, j0, nx, ny)
                                  : z[linear_index(max(i0 - 1, 0), j0, nx)];
  float z_avg = 0.25f * (z_north + z_south + z_east + z_west);
  float crest = max(0.f, z0 - z_avg);
  float eff_wind_speed = wind_speed * (1.0f + crest_speedup * crest * (float)nx);

  // current saltation trajectory state
  float current_x = (float)i0;
  float current_y = (float)j0;
  float current_z = z0;
  float angle = wind_angle;

  bool deposited = false;

  for (int hop = 0; hop < max_hops; ++hop)
  {
    // apply stochastic wind jitter
    float hop_angle = angle + (rand(&rng) * 2.f - 1.f) * jitter_angle;
    float2 wdir = (float2)(cos(hop_angle), sin(hop_angle));

    // base hop
    float base_hop = hop_length * eff_wind_speed;
    float next_x = current_x + base_hop * wdir.x;
    float next_y = current_y + base_hop * wdir.y;

    // handle domain wrapping / boundary
    if (periodic != 0)
    {
      next_x = fmod(fmod(next_x, (float)nx) + (float)nx, (float)nx);
      next_y = fmod(fmod(next_y, (float)ny) + (float)ny, (float)ny);
    }
    else
    {
      if (next_x < 0.f || next_x >= (float)nx || next_y < 0.f || next_y >= (float)ny)
      {
        // left domain: if inflow rate active, respawn at upwind edge, else lost
        if (sand_inflow_rate > 0.f && rand(&rng) < sand_inflow_rate)
        {
          next_x = (wdir.x >= 0.f) ? 0.f : (float)(nx - 1);
          next_y = rand(&rng) * (float)(ny - 1);
        }
        else
        {
          break;
        }
      }
    }

    int lx = clamp((int)floor(next_x + 0.5f), 0, nx - 1);
    int ly = clamp((int)floor(next_y + 0.5f), 0, ny - 1);
    int l_idx = linear_index(lx, ly, nx);
    float z_land = z[l_idx];

    // elevation-aware hop adjustment: downward gap carries particle further
    if (current_z > z_land)
    {
      float elev_drop = (current_z - z_land) * (float)nx;
      float extra_hop = elevation_hop_factor * elev_drop;
      next_x = current_x + (base_hop + extra_hop) * wdir.x;
      next_y = current_y + (base_hop + extra_hop) * wdir.y;

      if (periodic != 0)
      {
        next_x = fmod(fmod(next_x, (float)nx) + (float)nx, (float)nx);
        next_y = fmod(fmod(next_y, (float)ny) + (float)ny, (float)ny);
      }
      else if (next_x < 0.f || next_x >= (float)nx || next_y < 0.f || next_y >= (float)ny)
      {
        break;
      }

      lx = clamp((int)floor(next_x + 0.5f), 0, nx - 1);
      ly = clamp((int)floor(next_y + 0.5f), 0, ny - 1);
      l_idx = linear_index(lx, ly, nx);
      z_land = z[l_idx];
    }

    // shadow zone check using actual shadow talus
    bool in_shadow = helper_is_in_shadow(z, lx, ly, z_land, wdir, shadow_talus, nx, ny, periodic);

    bool should_deposit = false;
    if (in_shadow)
    {
      should_deposit = true;
    }
    else
    {
      bool has_sand = true;
      if (has_bedrock != 0)
      {
        has_sand = (z_land > bedrock[l_idx] + 1e-4f);
      }

      float p_dep = has_sand ? p_sand : p_bare;
      should_deposit = (rand(&rng) < p_dep);
    }

    if (should_deposit)
    {
      atomic_add_float(&z[l_idx], slab);
      deposited = true;
      break;
    }

    // continue hopping
    current_x = next_x;
    current_y = next_y;
    current_z = z_land;
  }

  // if particle didn't settle after max_hops, deposit at current position to conserve mass
  if (!deposited)
  {
    int lx = clamp((int)floor(current_x + 0.5f), 0, nx - 1);
    int ly = clamp((int)floor(current_y + 0.5f), 0, ny - 1);
    int l_idx = linear_index(lx, ly, nx);
    atomic_add_float(&z[l_idx], slab);
  }
}

void kernel sand_avalanche_step(global const float *z_in,
                                global float       *z_out,
                                global const float *bedrock,
                                const int           nx,
                                const int           ny,
                                const float         talus,
                                const float         collapse_rate,
                                const int           has_bedrock,
                                const int           periodic)
{
  int gx = get_global_id(0);
  int gy = get_global_id(1);

  if (gx >= nx || gy >= ny) return;

  // non-periodic boundaries: maintain closed boundary with extrapolation
  if (periodic == 0)
  {
    if (gx <= 0 || gx >= nx - 1 || gy <= 0 || gy >= ny - 1)
    {
      apply_boundaries_io(z_in, z_out, gx, gy, nx, ny);
      return;
    }
  }

  int idx_c = linear_index(gx, gy, nx);
  float z_c = z_in[idx_c];
  float sand_avail_c = (has_bedrock != 0) ? max(0.f, z_c - bedrock[idx_c]) : 1e10f;

  // 1. Calculate outgoing flux from (gx, gy) to lower neighbors exceeding talus
  float excess[8];
  float excess_sum = 0.f;
  float excess_max = 0.f;

  for (int k = 0; k < 8; ++k)
  {
    int ni = gx + d8_di[k];
    int nj = gy + d8_dj[k];

    float z_n;
    if (periodic != 0)
    {
      z_n = helper_sample_wrapped_height(z_in, ni, nj, nx, ny);
    }
    else
    {
      if (ni <= 0 || ni >= nx - 1 || nj <= 0 || nj >= ny - 1)
      {
        excess[k] = 0.f;
        continue;
      }
      z_n = z_in[linear_index(ni, nj, nx)];
    }

    float diff = z_c - z_n;
    float max_diff = talus * d8_dist[k];
    float e = diff - max_diff;

    if (e > 0.f)
    {
      excess[k] = e;
      excess_sum += e;
      excess_max = max(excess_max, e);
    }
    else
    {
      excess[k] = 0.f;
    }
  }

  float outflow = 0.f;
  if (excess_sum > 0.f && sand_avail_c > 0.f)
  {
    float desired_outflow = collapse_rate * 0.5f * excess_max;
    outflow = min(desired_outflow, sand_avail_c);
  }

  // 2. Calculate incoming flux from higher neighbors to (gx, gy)
  float inflow = 0.f;

  for (int k = 0; k < 8; ++k)
  {
    int ni = gx + d8_di[k];
    int nj = gy + d8_dj[k];

    int rx = ni;
    int ry = nj;
    if (periodic != 0)
    {
      rx = helper_wrap_coord(ni, nx);
      ry = helper_wrap_coord(nj, ny);
    }
    else
    {
      if (ni <= 0 || ni >= nx - 1 || nj <= 0 || nj >= ny - 1) continue;
    }

    int idx_n = linear_index(rx, ry, nx);
    float z_n = z_in[idx_n];
    float n_sand_avail = (has_bedrock != 0) ? max(0.f, z_n - bedrock[idx_n]) : 1e10f;
    if (n_sand_avail <= 0.f) continue;

    // recompute neighbor's excess toward all ITS active neighbors
    float n_excess_sum = 0.f;
    float n_excess_max = 0.f;
    float n_excess_to_c = 0.f;

    for (int p = 0; p < 8; ++p)
    {
      int qi = rx + d8_di[p];
      int qj = ry + d8_dj[p];

      float z_q;
      if (periodic != 0)
      {
        z_q = helper_sample_wrapped_height(z_in, qi, qj, nx, ny);
      }
      else
      {
        if (qi <= 0 || qi >= nx - 1 || qj <= 0 || qj >= ny - 1) continue;
        z_q = z_in[linear_index(qi, qj, nx)];
      }

      float e = z_n - z_q - talus * d8_dist[p];
      if (e > 0.f)
      {
        n_excess_sum += e;
        n_excess_max = max(n_excess_max, e);

        // direction from neighbor back to center is d8_rev[k]
        if (p == d8_rev[k])
        {
          n_excess_to_c = e;
        }
      }
    }

    if (n_excess_to_c > 0.f && n_excess_sum > 0.f)
    {
      float n_desired_outflow = collapse_rate * 0.5f * n_excess_max;
      float n_actual_outflow = min(n_desired_outflow, n_sand_avail);
      inflow += n_actual_outflow * (n_excess_to_c / n_excess_sum);
    }
  }

  float z_new = z_c - outflow + inflow;
  if (has_bedrock != 0)
  {
    z_new = max(z_new, bedrock[idx_c]);
  }

  z_out[idx_c] = z_new;
}
)""
