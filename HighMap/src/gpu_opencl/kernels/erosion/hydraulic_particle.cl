R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#define VELOCITY_INIT 1.f
#define VELOCITY_MIN 0.01f
#define VOLUME_MIN 0.001f
#define SLOPE_MIN 0.0001f
#define MAX_IT 5000

// --- HELPERS

float helper_bilinear_interp(const float f00,
                             const float f10,
                             const float f01,
                             const float f11,
                             const float u,
                             const float v)
{
  float a10 = f10 - f00;
  float a01 = f01 - f00;
  float a11 = f11 - f10 - f01 + f00;

  return f00 + a10 * u + a01 * v + a11 * u * v;
}

inline float helper_sample_height(global float *z,
                                  int           i,
                                  int           j,
                                  int           nx,
                                  int           ny,
                                  float         u,
                                  float         v)
{
  if (i == nx - 1) i--;
  if (j == ny - 1) j--;

  return helper_bilinear_interp(z[linear_index(i, j, nx)],
                                z[linear_index(i + 1, j, nx)],
                                z[linear_index(i, j + 1, nx)],
                                z[linear_index(i + 1, j + 1, nx)],
                                u,
                                v);
}

inline void helper_bilinear_deposition(global float *z,
                                       int           i,
                                       int           j,
                                       int           nx,
                                       int           ny,
                                       float         u,
                                       float         v,
                                       float         amount)
{
  if (i == nx - 1) i--;
  if (j == ny - 1) j--;

  float d1 = (1.f - u) * (1.f - v);
  float d2 = u * (1.f - v);
  float d3 = (1.f - u) * v;
  float d4 = u * v;

  z[linear_index(i, j, nx)] -= amount * d1;
  z[linear_index(i + 1, j, nx)] -= amount * d2;
  z[linear_index(i, j + 1, nx)] -= amount * d3;
  z[linear_index(i + 1, j + 1, nx)] -= amount * d4;
}

inline void helper_bank_collapse(global float *z,
                                 global float *bedrock,
                                 int           i,
                                 int           j,
                                 int           nx,
                                 int           ny,
                                 float         talus_slope,
                                 float         collapse_rate,
                                 int           has_bedrock)
{
  if (collapse_rate <= 0.f) return;
  if (!is_inside(i, j, nx, ny)) return;

  int   idx_c = linear_index(i, j, nx);
  float z_c = z[idx_c];

  int   di[8] = {-1, 1, 0, 0, -1, 1, -1, 1};
  int   dj[8] = {0, 0, -1, 1, -1, -1, 1, 1};
  float inv_dist[8] =
      {1.f, 1.f, 1.f, 1.f, 0.70710678f, 0.70710678f, 0.70710678f, 0.70710678f};

  for (int k = 0; k < 8; ++k)
  {
    int ni = i + di[k];
    int nj = j + dj[k];

    if (ni >= 0 && ni < nx && nj >= 0 && nj < ny)
    {
      int   idx_n = linear_index(ni, nj, nx);
      float z_n = z[idx_n];
      float slope = (z_n - z_c) * inv_dist[k];

      if (slope > talus_slope)
      {
        float excess = (slope - talus_slope) * collapse_rate;
        if (has_bedrock != 0)
        {
          float max_avail = max(0.f, z_n - bedrock[idx_n]);
          excess = min(excess, max_avail);
        }
        z[idx_n] -= excess * 0.5f;
        z[idx_c] += excess * 0.5f;
        z_c += excess * 0.5f;
      }
    }
  }
}

// --- MAIN KERNEL

void kernel hydraulic_particle(global float *z_in,
                               global float *bedrock,
                               global float *moisture_map,
                               global float *elevation_shift,
                               const int     nx,
                               const int     ny,
                               const int     nparticles,
                               const uint    seed,
                               const float   c_capacity,
                               const float   c_erosion,
                               const float   c_deposition,
                               const float   c_inertia,
                               const float   c_gravity,
                               const float   drag_rate,
                               const float   evap_rate,
                               const float   talus_slope,
                               const float   collapse_rate,
                               const int     has_bedrock,
                               const int     has_moisture_map,
                               const int     has_elevation_shift)
{
  float dt = 1.f;

  int id = get_global_id(0); // particle id
  if (id >= nparticles) return;

  uint   rng = pcg_hash(seed + id * 2u);
  float2 pos = {(nx - 2) * rand(&rng), (ny - 2) * rand(&rng)};
  int    i, j;
  float  u, v;

  update_interp_param(pos, &i, &j, &u, &v);

  float vel = VELOCITY_INIT;
  float volume = has_moisture_map > 0 ? moisture_map[linear_index(i, j, nx)]
                                      : 1.f;

  float2 dir = {0.f, 0.f};

  float s = 0.f;
  float evap_factor = 1.f - evap_rate;
  int   count = 0;

  while (volume > 1e-3f && count < MAX_IT)
  {
    count++; // for pathological cases...

    update_interp_param(pos, &i, &j, &u, &v);

    // stop if the particle reaches the domain limits
    if (!is_inside_gap(i, j, nx, ny, 1)) break;

    float2 gz;

    // compute continuous gradient at particle position
    {
      int idx00 = linear_index(i, j, nx);
      int idx10 = linear_index(i + 1, j, nx);
      int idx01 = linear_index(i, j + 1, nx);
      int idx11 = linear_index(i + 1, j + 1, nx);

      float z00 = z_in[idx00];
      float z10 = z_in[idx10];
      float z01 = z_in[idx01];
      float z11 = z_in[idx11];

      if (has_elevation_shift)
      {
        z00 += elevation_shift[idx00];
        z10 += elevation_shift[idx10];
        z01 += elevation_shift[idx01];
        z11 += elevation_shift[idx11];
      }

      gz.x = (1.f - v) * (z10 - z00) + v * (z11 - z01);
      gz.y = (1.f - u) * (z01 - z00) + u * (z11 - z10);
    }

    // ensure minimum slope
    float gz_norm = length(gz);

    // Option 4: Across flats (||grad z|| -> 0), preserve existing momentum dir
    // rather than resetting or artificially locking to a default vector
    float2 downhill_dir;
    if (gz_norm > SLOPE_MIN)
    {
      downhill_dir = (float2)(-gz.x / gz_norm, -gz.y / gz_norm);
    }
    else
    {
      downhill_dir = (length(dir) > 1e-4f) ? dir : (float2)(0.f, 0.f);
      gz_norm = SLOPE_MIN;
    }

    // Blend downhill direction with particle momentum
    // If terrain is flat, inertia dominates so particle carries smoothly across
    float eff_inertia = (gz_norm <= 2.f * SLOPE_MIN) ? max(c_inertia, 0.9f)
                                                     : c_inertia;
    dir = mix(downhill_dir, dir, eff_inertia);

    if (length(dir) > 0.f)
      dir /= length(dir);
    else
      dir = (float2)(0.f, 0.f);

    if (vel < VELOCITY_MIN) break;

    float zp = helper_sample_height(z_in, i, j, nx, ny, u, v);

    // backup previous position
    int   ip = i;
    int   jp = j;
    float up = u;
    float vp = v;

    // move particle
    pos += dt * dir;

    // elevation at new position
    update_interp_param(pos, &i, &j, &u, &v);
    if (!is_inside(i, j, nx, ny)) break;

    float z = helper_sample_height(z_in, i, j, nx, ny, u, v);
    float dz = zp - z;
    // Elevation drop and stream power capacity
    float sc = max(0.f, c_capacity * volume * vel * dz);
    float delta_sc = dt * (sc - s);
    float amount = 0.f;

    // Sample local sediment thickness if bedrock is present
    float local_sed = 0.f;
    if (has_bedrock != 0)
    {
      float zb = helper_sample_height(bedrock, i, j, nx, ny, u, v);
      local_sed = max(0.f, z - zb);
    }

    // Deposition / Erosion
    if (delta_sc < 0.f || dz < 0.f)
    {
      // Deposition: amount is negative so helper_bilinear_deposition adds
      // height to z_in, and adding amount to s decreases sediment carried.
      float to_deposit = (dz < 0.f) ? max(-dz, -delta_sc) : -delta_sc;
      amount = -c_deposition * min(s, to_deposit);
      helper_bilinear_deposition(z_in, ip, jp, nx, ny, up, vp, amount);
    }
    else
    {
      // Erosion: amount is positive so helper_bilinear_deposition subtracts
      // height from z_in.
      // Soft sediment / deposited alluvium has higher erodibility than hard
      // bedrock.
      float eff_c_erosion = c_erosion;
      if (has_bedrock != 0 && local_sed > 1e-4f)
      {
        // Boost erosion rate on loose deposited sediment so depressions are
        // easily flushed
        eff_c_erosion = min(1.f, c_erosion * 4.f);
      }

      float max_erode = (dz > 0.f) ? dz : delta_sc;
      amount = eff_c_erosion * min(delta_sc, max_erode);
      helper_bilinear_deposition(z_in, ip, jp, nx, ny, up, vp, amount);
    }

    s += amount;

    // Bedrock limit enforcement across interpolated footprint
    if (amount > 0.f && has_bedrock != 0)
    {
      int i0 = ip, j0 = jp;
      if (i0 == nx - 1) i0--;
      if (j0 == ny - 1) j0--;
      int idx00 = linear_index(i0, j0, nx);
      int idx10 = linear_index(i0 + 1, j0, nx);
      int idx01 = linear_index(i0, j0 + 1, nx);
      int idx11 = linear_index(i0 + 1, j0 + 1, nx);
      z_in[idx00] = max(bedrock[idx00], z_in[idx00]);
      z_in[idx10] = max(bedrock[idx10], z_in[idx10]);
      z_in[idx01] = max(bedrock[idx01], z_in[idx01]);
      z_in[idx11] = max(bedrock[idx11], z_in[idx11]);
    }

    // bank collapse across the bilinear footprint
    if (amount > 0.f && collapse_rate > 0.f)
    {
      int i0 = ip;
      int j0 = jp;
      if (i0 == nx - 1) i0--;
      if (j0 == ny - 1) j0--;

      float d1 = (1.f - up) * (1.f - vp);
      float d2 = up * (1.f - vp);
      float d3 = (1.f - up) * vp;
      float d4 = up * vp;

      if (d1 > 0.f)
        helper_bank_collapse(z_in,
                             bedrock,
                             i0,
                             j0,
                             nx,
                             ny,
                             talus_slope,
                             collapse_rate * d1,
                             has_bedrock);
      if (d2 > 0.f)
        helper_bank_collapse(z_in,
                             bedrock,
                             i0 + 1,
                             j0,
                             nx,
                             ny,
                             talus_slope,
                             collapse_rate * d2,
                             has_bedrock);
      if (d3 > 0.f)
        helper_bank_collapse(z_in,
                             bedrock,
                             i0,
                             j0 + 1,
                             nx,
                             ny,
                             talus_slope,
                             collapse_rate * d3,
                             has_bedrock);
      if (d4 > 0.f)
        helper_bank_collapse(z_in,
                             bedrock,
                             i0 + 1,
                             j0 + 1,
                             nx,
                             ny,
                             talus_slope,
                             collapse_rate * d4,
                             has_bedrock);
    }

    // Velocity update:
    // Option 4: Soften uphill deceleration for low-angle depressions/obstacles
    // (scale uphill gravity loss with local adverse slope ratio)
    float g_eff = c_gravity;
    if (dz < 0.f)
    {
      // Shallow uphill rises (e.g. -dz < 0.05) face significantly reduced
      // deceleration, allowing incoming fast water to glide over low-angle
      // sills.
      float slope_atten = clamp(-dz * 20.f, 0.2f, 1.f);
      g_eff = c_gravity * slope_atten;
    }
    vel = sqrt(max(0.f, vel * vel - dz * g_eff));

    // Damping drag
    vel *= (1.f - drag_rate);

    volume *= evap_factor;
  }
}
)""
