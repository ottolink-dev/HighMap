R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel mountain_range_radial(global float *output,
                                  global float *ctrl_param,
                                  global float *noise_x,
                                  global float *noise_y,
                                  global float *angle,
                                  const int     nx,
                                  const int     ny,
                                  const float   kx,
                                  const float   ky,
                                  const uint    seed,
                                  const float   half_width,
                                  const float   angle_spread_ratio,
                                  const float   core_size_ratio,
                                  const float2  center,
                                  const int     octaves,
                                  const float   weight,
                                  const float   persistence,
                                  const float   lacunarity,
                                  const int     has_ctrl_param,
                                  const int     has_noise_x,
                                  const int     has_noise_y,
                                  const int     has_angle,
                                  const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float ct = has_ctrl_param > 0 ? ctrl_param[index] : 1.f;
  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  // overall amplitude (Gaussian pulse)
  float2 xy = g_to_xy(g, nx, ny, 1.f, 1.f, 0.f, 0.f, bbox) - center;
  float  r2 = dot(xy, xy);
  float  amp = exp(-0.5f * r2 / (half_width * half_width));

  // noise angle perpendicular to radius
  float  theta = atan2(xy.y, xy.x) + 1.57080f;
  float2 dir = {cos(theta), sin(theta)};

  // align roughness with amplitude
  ct *= amp;
  float noise = gabor_wave_scalar_fbm(pos,
                                      dir,
                                      angle_spread_ratio,
                                      octaves,
                                      (1.f - ct) + ct * weight,
                                      persistence,
                                      lacunarity,
                                      fseed);

  // smoothing at origin to avoid numerical artifacts
  float r2_max = core_size_ratio / max(kx, ky);
  float t = min(1.f, r2 / r2_max);
  t = sqrt(t) * (1.f - exp(-500.f * t));
  t = smoothstep3_upper(t);

  output[index] = amp * lerp(1.f, 0.5f * noise + 0.5f, t);

  if (has_angle) angle[index] = theta;
}
)""
