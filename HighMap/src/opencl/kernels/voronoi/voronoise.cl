R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// Adapted from https://www.shadertoy.com/view/Xd23Dh
// Original author: Inigo Quilez
// License: MIT (see THIRD_PARTY_LICENSES.md)

float3 helper_voronoise_hash3(const float2 p, float fseed)
{
  float3 q = (float3)(dot(p, (float2)(127.1f, 311.7f)),
                      dot(p, (float2)(269.5f, 183.3f)),
                      dot(p, (float2)(419.2f, 371.9f)));
  float3 qi;
  return fract(sin(q) * 43758.5453f + fseed, &qi);
}

float noise_voronoise(const float2 p,
                      const float  u_param,
                      const float  v_param,
                      const float  fseed)
{
  float k = 1.f + 63.f * pow_float(1.f - v_param, 6.f);

  float2 i = floor(p);
  float2 pi;
  float2 f = fract(p, &pi);

  float2 a = (float2)(0.f, 0.f);
  for (int y = -2; y <= 2; y++)
    for (int x = -2; x <= 2; x++)
    {
      float2 g = (float2)(x, y);
      float3 o = helper_voronoise_hash3(i + g, fseed) *
                 (float3)(u_param, u_param, 1.f);
      float2 d = g - f + o.xy;
      float  w = pow_float(1.f - smoothstep(0.f, 1.414f, length(d)), k);
      a += (float2)(o.z * w, w);
    }

  return a.x / a.y;
}
// <<<

float noise_voronoise_fbm(const float2 p,
                          const float  u_param,
                          const float  v_param,
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
    float v = noise_voronoise(p * nf, u_param, v_param, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(2.f * v, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

void kernel voronoise(global float *output,
                      global float *noise_x,
                      global float *noise_y,
                      const int     nx,
                      const int     ny,
                      const float   kx,
                      const float   ky,
                      const float   u_param,
                      const float   v_param,
                      const uint    seed,
                      const int     has_noise_x,
                      const int     has_noise_y,
                      const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  output[index] = noise_voronoise(pos, u_param, v_param, fseed);
}

void kernel voronoise_fbm(global float *output,
                          global float *ctrl_param,
                          global float *noise_x,
                          global float *noise_y,
                          const int     nx,
                          const int     ny,
                          const float   kx,
                          const float   ky,
                          const float   u_param,
                          const float   v_param,
                          const uint    seed,
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

  output[index] = noise_voronoise_fbm(pos,
                                      u_param,
                                      v_param,
                                      octaves,
                                      (1.f - ct) + ct * weight,
                                      persistence,
                                      lacunarity,
                                      fseed);
}
)""
