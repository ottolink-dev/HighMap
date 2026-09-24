R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

float sample_angle_field(global const float *angle,
                         const float2        cell_pos,
                         const int           nx,
                         const int           ny,
                         const float         kx,
                         const float         ky,
                         const float4        bbox)
{
  float u = (cell_pos.x / kx - bbox.x) / (bbox.y - bbox.x);
  float v = (cell_pos.y / ky - bbox.z) / (bbox.w - bbox.z);

  float px = clamp(u * (float)nx - 0.5f, 0.f, (float)(nx - 1));
  float py = clamp(v * (float)ny - 0.5f, 0.f, (float)(ny - 1));

  int x0 = (int)px;
  int y0 = (int)py;
  int x1 = min(x0 + 1, nx - 1);
  int y1 = min(y0 + 1, ny - 1);

  float fx = px - (float)x0;
  float fy = py - (float)y0;

  float a00 = angle[linear_index(x0, y0, nx)];
  float a10 = angle[linear_index(x1, y0, nx)];
  float a01 = angle[linear_index(x0, y1, nx)];
  float a11 = angle[linear_index(x1, y1, nx)];

  float w00 = (1.f - fx) * (1.f - fy);
  float w10 = fx * (1.f - fy);
  float w01 = (1.f - fx) * fy;
  float w11 = fx * fy;

  float cos_val = w00 * cos(a00) + w10 * cos(a10) + w01 * cos(a01) +
                  w11 * cos(a11);
  float sin_val = w00 * sin(a00) + w10 * sin(a10) + w01 * sin(a01) +
                  w11 * sin(a11);

  return atan2(sin_val, cos_val);
}

float2 base_phase_field(global const float *angle_field,
                        const float2        p,
                        const float2        jitter,
                        const float         freq,
                        const float         normalization,
                        const float         fseed,
                        const int           nx,
                        const int           ny,
                        const float         kx,
                        const float         ky,
                        const float4        bbox)
{
  float2 p_int = floor(p);
  float2 p_frac = p - p_int;

  float2      phase_dir = 0.f;
  float       weight_sum = 0.f;
  const float tau = 6.283185307179586f;

  for (int i = -1; i <= 2; ++i)
  {
    for (int j = -1; j <= 2; ++j)
    {
      float2 grid_offset = (float2)(i, j);
      float2 grid_point = p_int + grid_offset;

      // random offset for the cell point between -0.5 and 0.5 on each axis
      float2 random_offset = (-1.f + 2.f * hash22f_poly(grid_point, fseed)) *
                             0.5f * jitter;
      float2 cell_point = grid_point + random_offset;

      float cell_angle =
          sample_angle_field(angle_field, cell_point, nx, ny, kx, ky, bbox);

      // wave vector along the guidance direction
      float2 wave_vec = (float2)(cos(cell_angle), sin(cell_angle)) *
                        (freq * tau);

      float2 vec_from_cell = p_frac - grid_offset - random_offset;

      // bell-shaped weight function: max(0, exp(-2 * r^2) - 0.01111)
      float sqr_dist = dot(vec_from_cell, vec_from_cell);
      float weight = exp(-sqr_dist * 2.f);
      weight = max(0.f, weight - 0.01111f);

      weight_sum += weight;

      float wave_input = dot(vec_from_cell, wave_vec);

      phase_dir += (float2)(cos(wave_input), sin(wave_input)) * weight;
    }
  }

  if (weight_sum <= 1e-6f)
  {
    return 0.f;
  }

  float2 interpolated = phase_dir / weight_sum;
  float  magnitude = length(interpolated);
  float  clamped_mag = max(1.f - normalization, magnitude);

  if (clamped_mag > 1e-6f)
  {
    return interpolated / clamped_mag;
  }
  return interpolated;
}

void kernel phase_field(global float *angle,
                        global float *phase,
                        global float *ctrl_param,
                        global float *noise_x,
                        global float *noise_y,
                        global float *modulus,
                        const int     nx,
                        const int     ny,
                        const float   kx,
                        const float   ky,
                        const uint    seed,
                        const float2  jitter,
                        const float   normalization,
                        const float   kp,
                        const int     has_ctrl_param,
                        const int     has_noise_x,
                        const int     has_noise_y,
                        const int     has_modulus,
                        const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int idx = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float ct = has_ctrl_param > 0 ? ctrl_param[idx] : 1.f;
  float dx = has_noise_x > 0 ? noise_x[idx] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[idx] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  float2 phacelle = base_phase_field(angle,
                                     pos,
                                     jitter,
                                     ct * kp,
                                     normalization,
                                     fseed,
                                     nx,
                                     ny,
                                     kx,
                                     ky,
                                     bbox);

  phase[idx] = atan2(phacelle.y, phacelle.x);

  if (has_modulus) modulus[idx] = length(phacelle);
}
)""
