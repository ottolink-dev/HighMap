R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

float wavelet_noise_base(const float2 p,
                         const float  phase,
                         const float  kw_multiplier,
                         const float  vorticity,
                         const float  density,
                         const int    octaves,
                         const float  weight,
                         const float  persistence,
                         const float  lacunarity,
                         const float  fseed)
{
  // Adapted from https://www.shadertoy.com/view/wsBfzK
  // Original author:  Martijn Steinrucken
  // License: MIT (see THIRD_PARTY_LICENSES.md)

  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;

  float  a = 0.0f;
  float2 pm = p;

  for (int it = 0; it < octaves; it++)
  {
    float2 q = pm * nf;

    float qrd = hash12f(floor(q), fseed);

    a = qrd * 1e3;
    a += phase * qrd * vorticity;

    // q = (fract(q) - 0.5) * rotate2d(a);
    float2 qi;
    float2 fq = fract(q, &qi) - (float2)(0.5f, 0.5f);
    q = rotate2d(fq, a);

    float dq = dot(q, q);

    float v = (0.5f * sin(6.2832f * q.x * kw_multiplier + phase) + 0.5f) *
              smoothstep3_gl(dq, 0.25f, 0.0f);

    // randomly choose to not fill a cell (based on the density)
    if (qrd > density) v = 0.f;

    // pm = pm * mat2(0.54,-0.84, 0.84, 0.54) + i;
    float2 rot;
    rot.x = pm.x * 0.54f + pm.y * 0.84f;
    rot.y = -pm.x * 0.84f + pm.y * 0.54f;
    pm = rot + (float)it;

    // std fbm
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }

  return n;
}

void kernel wavelet_noise(global float *output,
                          global float *ctrl_param,
                          global float *noise_x,
                          global float *noise_y,
                          const int     nx,
                          const int     ny,
                          const float   kx,
                          const float   ky,
                          const uint    seed,
                          const float   kw_multiplier,
                          const float   vorticity,
                          const float   density,
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

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  output[index] = wavelet_noise_base(pos,
                                     ct,
                                     kw_multiplier,
                                     vorticity,
                                     density,
                                     octaves,
                                     weight,
                                     persistence,
                                     lacunarity,
                                     fseed);
}
)""
