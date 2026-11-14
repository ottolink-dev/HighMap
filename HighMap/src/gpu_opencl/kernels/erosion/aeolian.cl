R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

/* Clean-room GPU implementation of the aeolian desert simulation algorithm
 * based on Paris et al. (Pacific Graphics 2019 / CGF), "Desertscapes
 * Simulation". Models wind-driven sand transport (saltation & reptation), wind
 * shadowing, bedrock abrasion (yardang formation), vegetation stabilization
 * (nabkhas), and sediment angle-of-repose relaxation. */

#define AEOLIAN_MAX_BOUNCE 3

// --- Helpers

inline float helper_aeolian_height(global const float *bedrock,
                                   global const float *sediment,
                                   int                 i,
                                   int                 j,
                                   int                 nx,
                                   int                 ny)
{
  int idx = linear_index(clamp(i, 0, nx - 1), clamp(j, 0, ny - 1), nx);
  return bedrock[idx] + sediment[idx];
}

inline float helper_aeolian_height_bilinear(global const float *bedrock,
                                            global const float *sediment,
                                            float               x,
                                            float               y,
                                            int                 nx,
                                            int                 ny)
{
  // Toroidal coordinates
  float tx = fmod(fmod(x, (float)nx) + (float)nx, (float)nx);
  float ty = fmod(fmod(y, (float)ny) + (float)ny, (float)ny);

  int x0 = (int)floor(tx);
  int y0 = (int)floor(ty);
  int x1 = (x0 + 1) % nx;
  int y1 = (y0 + 1) % ny;

  float fx = tx - (float)x0;
  float fy = ty - (float)y0;

  int idx00 = linear_index(x0, y0, nx);
  int idx10 = linear_index(x1, y0, nx);
  int idx01 = linear_index(x0, y1, nx);
  int idx11 = linear_index(x1, y1, nx);

  float h00 = bedrock[idx00] + sediment[idx00];
  float h10 = bedrock[idx10] + sediment[idx10];
  float h01 = bedrock[idx01] + sediment[idx01];
  float h11 = bedrock[idx11] + sediment[idx11];

  float h0 = mix(h00, h10, fx);
  float h1 = mix(h01, h11, fx);

  return mix(h0, h1, fy);
}

inline float helper_aeolian_is_in_shadow(global const float *bedrock,
                                         global const float *sediment,
                                         float               curr_x,
                                         float               curr_y,
                                         int                 nx,
                                         int                 ny,
                                         float2              wind_dir,
                                         float               cell_size,
                                         float               tan_shadow_min,
                                         float               tan_shadow_max,
                                         float               max_shadow_dist)
{
  if (length(wind_dir) < 1e-5f) return 0.f;

  float2 uwind = normalize(wind_dir);
  float  curr_h = helper_aeolian_height_bilinear(bedrock,
                                                sediment,
                                                curr_x,
                                                curr_y,
                                                nx,
                                                ny);

  float max_shadow = 0.f;

  for (float step = 1.f; step <= max_shadow_dist; step += 1.f)
  {
    // raymarch upwind along continuous direction
    float up_x = curr_x - step * uwind.x;
    float up_y = curr_y - step * uwind.y;

    float up_h = helper_aeolian_height_bilinear(bedrock,
                                                sediment,
                                                up_x,
                                                up_y,
                                                nx,
                                                ny);
    float dh = up_h - curr_h;
    float dist_phys = step * cell_size;
    float tan_angle = dh / dist_phys;

    float shadow = smoothstep(tan_shadow_min, tan_shadow_max, tan_angle);
    max_shadow = max(max_shadow, shadow);
    if (max_shadow >= 1.f) break;
  }

  return max_shadow;
}

// Bilinear lookup of a precomputed shadow map (values sampled at cell centers)
inline float helper_aeolian_shadow_lookup(global const float *shadow_map,
                                          float               x,
                                          float               y,
                                          int                 nx,
                                          int                 ny)
{
  float tx = fmod(fmod(x - 0.5f, (float)nx) + (float)nx, (float)nx);
  float ty = fmod(fmod(y - 0.5f, (float)ny) + (float)ny, (float)ny);

  int x0 = (int)floor(tx);
  int y0 = (int)floor(ty);
  int x1 = (x0 + 1) % nx;
  int y1 = (y0 + 1) % ny;

  float fx = tx - (float)x0;
  float fy = ty - (float)y0;

  float s0 = mix(shadow_map[linear_index(x0, y0, nx)],
                 shadow_map[linear_index(x1, y0, nx)],
                 fx);
  float s1 = mix(shadow_map[linear_index(x0, y1, nx)],
                 shadow_map[linear_index(x1, y1, nx)],
                 fx);

  return mix(s0, s1, fy);
}

// --- Shadow map kernel (one raymarch per cell instead of one per walker hop)

void kernel aeolian_shadow_map(global const float *bedrock,
                               global const float *sediment,
                               global float       *shadow_map,
                               const int           nx,
                               const int           ny,
                               const float         wind_x,
                               const float         wind_y,
                               const float         cell_size,
                               const float         tan_shadow_min,
                               const float         tan_shadow_max,
                               const float         max_shadow_dist)
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  if (i >= nx || j >= ny) return;

  shadow_map[linear_index(i, j, nx)] = helper_aeolian_is_in_shadow(
      bedrock,
      sediment,
      (float)i + 0.5f,
      (float)j + 0.5f,
      nx,
      ny,
      (float2)(wind_x, wind_y),
      cell_size,
      tan_shadow_min,
      tan_shadow_max,
      max_shadow_dist);
}

// --- Main saltation & abrasion kernel

void kernel aeolian_transport(global float       *bedrock,
                              global float       *sediment,
                              global const float *vegetation,
                              global const float *abrasion_mask,
                              const int           nx,
                              const int           ny,
                              const int           n_samples,
                              const uint          seed,
                              const uint          step,
                              const float         wind_x,
                              const float         wind_y,
                              const float         matter_to_move,
                              const float         cell_size,
                              const float         tan_shadow_min,
                              const float         tan_shadow_max,
                              const float         max_shadow_dist,
                              const float         abrasion_rate,
                              const int           has_vegetation,
                              const int           has_abrasion_mask,
                              const int           use_shadow_map,
                              global const float *shadow_map)
{
  int id = get_global_id(0);
  if (id >= n_samples) return;

  uint rng = pcg_hash(seed + id * 4u + step * 0x9e3779b9u);

  // 1. Pick a random grid position
  int start_i = (int)(rand(&rng) * (float)nx);
  int start_j = (int)(rand(&rng) * (float)ny);
  start_i = clamp(start_i, 0, nx - 1);
  start_j = clamp(start_j, 0, ny - 1);

  int start_idx = linear_index(start_i, start_j, nx);

  // Check if sediment is available to lift
  if (sediment[start_idx] < matter_to_move) return;

  float2 base_wind = (float2)(wind_x, wind_y);
  float  base_wind_len = length(base_wind);

  // Modulate wind by local terrain gradient projected onto wind direction
  float h_east = helper_aeolian_height(bedrock,
                                       sediment,
                                       (start_i + 1) % nx,
                                       start_j,
                                       nx,
                                       ny);
  float h_west = helper_aeolian_height(bedrock,
                                       sediment,
                                       (start_i - 1 + nx) % nx,
                                       start_j,
                                       nx,
                                       ny);
  float h_north = helper_aeolian_height(bedrock,
                                        sediment,
                                        start_i,
                                        (start_j + 1) % ny,
                                        nx,
                                        ny);
  float h_south = helper_aeolian_height(bedrock,
                                        sediment,
                                        start_i,
                                        (start_j - 1 + ny) % ny,
                                        nx,
                                        ny);

  float  inv_2dx = 0.5f / cell_size;
  float2 grad = (float2)((h_east - h_west) * inv_2dx,
                         (h_north - h_south) * inv_2dx);
  float  wind_speedup = 1.f;

  if (base_wind_len > 1e-5f)
  {
    float2 u_wind = base_wind / base_wind_len;
    float  slope_along_wind = dot(grad, u_wind);
    wind_speedup = 1.f + 0.5f * clamp(slope_along_wind, 0.f, 1.f);
  }

  float2 eff_wind = base_wind * wind_speedup;

  // Wind shadow check at lifting site: sheltered cells resist lifting
  float2 curr_pos = (float2)((float)start_i + 0.5f, (float)start_j + 0.5f);
  float  shadow_prob = 0.f;

  if (use_shadow_map != 0)
    shadow_prob = helper_aeolian_shadow_lookup(shadow_map,
                                               curr_pos.x,
                                               curr_pos.y,
                                               nx,
                                               ny);
  else
    shadow_prob = helper_aeolian_is_in_shadow(bedrock,
                                              sediment,
                                              curr_pos.x,
                                              curr_pos.y,
                                              nx,
                                              ny,
                                              eff_wind,
                                              cell_size,
                                              tan_shadow_min,
                                              tan_shadow_max,
                                              max_shadow_dist);
  if (rand(&rng) < shadow_prob) return;

  // Vegetation protection at lifting site
  if (has_vegetation != 0)
  {
    float veg_val = vegetation[start_idx];
    if (rand(&rng) < veg_val) return;
  }

  // 2. Saltation hop loop in continuous subpixel coordinates
  int bounce = 0;
  int dest_i = start_i;
  int dest_j = start_j;

  while (bounce < AEOLIAN_MAX_BOUNCE)
  {
    // Hop downwind along continuous wind vector
    curr_pos += eff_wind;

    // Toroidal wrap boundaries
    curr_pos.x = fmod(fmod(curr_pos.x, (float)nx) + (float)nx, (float)nx);
    curr_pos.y = fmod(fmod(curr_pos.y, (float)ny) + (float)ny, (float)ny);

    dest_i = clamp((int)floor(curr_pos.x), 0, nx - 1);
    dest_j = clamp((int)floor(curr_pos.y), 0, ny - 1);

    int dest_idx = linear_index(dest_i, dest_j, nx);

    // Bedrock abrasion: happens when saltating sand impacts exposed bedrock
    if (abrasion_rate > 0.f && sediment[dest_idx] < 0.5f * matter_to_move)
    {
      if (rand(&rng) < 0.25f)
      {
        float rock_susceptibility = 1.f;
        if (has_abrasion_mask != 0)
          rock_susceptibility = abrasion_mask[dest_idx];

        float veg_prot = (has_vegetation != 0) ? vegetation[dest_idx] : 0.f;
        float abr = abrasion_rate * rock_susceptibility * (1.f - veg_prot) *
                    length(eff_wind);
        if (abr > 0.f) bedrock[dest_idx] -= abr;
      }
    }

    // Deposition probability test
    float p = rand(&rng);
    float dest_shadow = 0.f;

    if (use_shadow_map != 0)
      dest_shadow = helper_aeolian_shadow_lookup(shadow_map,
                                                 curr_pos.x,
                                                 curr_pos.y,
                                                 nx,
                                                 ny);
    else
      dest_shadow = helper_aeolian_is_in_shadow(bedrock,
                                                sediment,
                                                curr_pos.x,
                                                curr_pos.y,
                                                nx,
                                                ny,
                                                eff_wind,
                                                cell_size,
                                                tan_shadow_min,
                                                tan_shadow_max,
                                                max_shadow_dist);

    float veg_dest = (has_vegetation != 0) ? vegetation[dest_idx] : 0.f;
    float dep_threshold = (sediment[dest_idx] > 0.f) ? (0.6f + 0.4f * veg_dest)
                                                     : (0.4f + 0.6f * veg_dest);

    if (p < dest_shadow || p < dep_threshold)
    {
      break;
    }

    bounce++;
  }

  int dest_idx = linear_index(dest_i, dest_j, nx);

  if (dest_idx != start_idx)
  {
    sediment[start_idx] -= matter_to_move;
    sediment[dest_idx] += matter_to_move;
  }
}

// --- Angle of repose stabilization (double-buffered conservative relaxation)

void kernel aeolian_stabilize_sediment(global const float *bedrock,
                                       global const float *sed_in,
                                       global float       *sed_out,
                                       const int           nx,
                                       const int           ny,
                                       const float         cell_size,
                                       const float         tan_repose,
                                       const float         rate)
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  if (i >= nx || j >= ny) return;

  int   idx = linear_index(i, j, nx);
  float s_c = sed_in[idx];
  float b_c = bedrock[idx];
  float h_c = b_c + s_c;

  if (s_c <= 0.f)
  {
    sed_out[idx] = s_c;
    return;
  }

  // D8 neighbor offsets and physical distances
  int   di[8] = {-1, 1, 0, 0, -1, 1, -1, 1};
  int   dj[8] = {0, 0, -1, 1, -1, -1, 1, 1};
  float dist[8] = {cell_size,
                   cell_size,
                   cell_size,
                   cell_size,
                   1.41421356f * cell_size,
                   1.41421356f * cell_size,
                   1.41421356f * cell_size,
                   1.41421356f * cell_size};

  // 1. Calculate outflow to neighbors that exceed angle of repose
  float excess[8];
  float excess_sum = 0.f;
  float excess_max = 0.f;

  for (int k = 0; k < 8; ++k)
  {
    int ni = (i + di[k] % nx + nx) % nx;
    int nj = (j + dj[k] % ny + ny) % ny;

    int   n_idx = linear_index(ni, nj, nx);
    float h_n = bedrock[n_idx] + sed_in[n_idx];
    float dh = h_c - h_n;
    float slope = dh / dist[k];

    if (slope > tan_repose)
    {
      float e = (slope - tan_repose) * dist[k];
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
  if (excess_sum > 0.f)
  {
    outflow = min(s_c, rate * 0.5f * excess_max);
  }

  // 2. Calculate inflow from neighbors
  float inflow = 0.f;

  for (int k = 0; k < 8; ++k)
  {
    int ni = (i + di[k] % nx + nx) % nx;
    int nj = (j + dj[k] % ny + ny) % ny;

    int   n_idx = linear_index(ni, nj, nx);
    float s_n = sed_in[n_idx];
    if (s_n <= 0.f) continue;

    float h_n = bedrock[n_idx] + s_n;

    // Check neighbor's slope toward center
    float dh_to_c = h_n - h_c;
    float slope_to_c = dh_to_c / dist[k];

    if (slope_to_c > tan_repose)
    {
      // Recompute neighbor's total excess to apportion fairly
      float n_excess_sum = 0.f;
      float n_excess_max = 0.f;
      float e_to_c = (slope_to_c - tan_repose) * dist[k];

      for (int p = 0; p < 8; ++p)
      {
        int   qi = (ni + di[p] % nx + nx) % nx;
        int   qj = (nj + dj[p] % ny + ny) % ny;
        int   q_idx = linear_index(qi, qj, nx);
        float h_q = bedrock[q_idx] + sed_in[q_idx];
        float n_slope = (h_n - h_q) / dist[p];

        if (n_slope > tan_repose)
        {
          float ne = (n_slope - tan_repose) * dist[p];
          n_excess_sum += ne;
          n_excess_max = max(n_excess_max, ne);
        }
      }

      if (n_excess_sum > 0.f)
      {
        float n_outflow = min(s_n, rate * 0.5f * n_excess_max);
        inflow += n_outflow * (e_to_c / n_excess_sum);
      }
    }
  }

  sed_out[idx] = max(0.f, s_c - outflow + inflow);
}
)""
