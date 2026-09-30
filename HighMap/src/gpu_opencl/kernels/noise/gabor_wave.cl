R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// Adapted from https://www.shadertoy.com/view/clGyWm
// Original author:  Inigo Quilez
// License: MIT (see THIRD_PARTY_LICENSES.md)

float gabor_wave_scalar(const float2 p,
                        const float2 dir,
                        const float  angle_spread_ratio,
                        const float  fseed)
{
  float2 ip = floor(p);
  float2 qi;
  float2 fp = fract(p, &qi);

  const float fr = 6.283185f;
  const float fa = 4.f;

  float av = 0.f;
  float at = 0.f;

  for (int j = -2; j <= 2; j++)
    for (int i = -2; i <= 2; i++)
    {
      float2 o = {i, j};
      float2 h = hash22f_poly(ip + o, fseed);
      float2 r = fp - (o + h);

      float2 s = {11.f, 31.f};
      float2 k = normalize(dir +
                           angle_spread_ratio *
                               (2.f * hash22f_poly(ip + o + s, fseed) - 1.f));

      float d = dot(r, r);
      float l = dot(r, k);
      float w = exp(-fa * d);
      float cs = cos(fr * l);

      av += w * cs;
      at += w;
    }

  return av / at;
}

// <<<

float gabor_wave_scalar_fbm(const float2 p,
                            const float2 dir,
                            const float  angle_spread_ratio,
                            const int    octaves,
                            const float  weight,
                            const float  persistence,
                            const float  lacunarity,
                            const float  fseed)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  for (int i = 0; i < octaves; i++)
  {
    float v = gabor_wave_scalar(p * nf, dir, angle_spread_ratio, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

void kernel gabor_wave(global float *output,
                       global float *angle,
                       const int     nx,
                       const int     ny,
                       const float   kx,
                       const float   ky,
                       const uint    seed,
                       const float   angle_spread_ratio,
                       const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  // "0.5f * kx" to keep it coherent with Perlin
  float2 pos = g_to_xy(g, nx, ny, 0.5f * kx, 0.5f * ky, 0.f, 0.f, bbox);

  float2 dir = angle_to_dir(angle[index]);

  output[index] = gabor_wave_scalar(pos, dir, angle_spread_ratio, fseed);
}

void kernel gabor_wave_fbm(global float *output,
                           global float *angle,
                           global float *ctrl_param,
                           global float *noise_x,
                           global float *noise_y,
                           const int     nx,
                           const int     ny,
                           const float   kx,
                           const float   ky,
                           const uint    seed,
                           const float   angle_spread_ratio,
                           const int     octaves,
                           const float   weight,
                           const float   persistence,
                           const float   lacunarity,
                           const int     has_ctrl_param,
                           const int     has_noise_x,
                           const int     has_noise_y,
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

  // "0.5f * kx" to keep it coherent with Perlin
  float2 pos = g_to_xy(g, nx, ny, 0.5f * kx, 0.5f * ky, dx, dy, bbox);

  float2 dir = angle_to_dir(angle[index]);

  output[index] = gabor_wave_scalar_fbm(pos,
                                        dir,
                                        angle_spread_ratio,
                                        octaves,
                                        (1.f - ct) + ct * weight,
                                        persistence,
                                        lacunarity,
                                        fseed);
}
)""
