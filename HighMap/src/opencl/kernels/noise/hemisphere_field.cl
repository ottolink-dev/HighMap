R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

float hf_base(const float2 pos,
              const float  rmin,
              const float  rmax,
              const float  amplitude_random_ratio,
              const float  density,
              const float2 jitter,
              const float  shift,
              const float  fseed)
{
  float  val = FLT_MAX;
  float2 pos_i = floor(pos);
  float  amp_rnd = 0.f;

  for (int dx = -2; dx <= 2; dx++)
    for (int dy = -2; dy <= 2; dy++)
    {
      float2 pos_nbrs = pos_i + (float2)(dx, dy);

      // occurence probability
      float does_spawn_proba = hash12f(pos_nbrs, fseed);

      if (does_spawn_proba < density)
      {
        float  radius = rmin + (rmax - rmin) * hash12f(pos_nbrs, fseed);
        float2 delta = jitter * hash22f(pos_nbrs, fseed);
        float2 center = pos_nbrs + delta;

        float r = length(pos - center) / radius;

        if (r < val)
        {
          val = r;
          amp_rnd = hash12f(pos_nbrs + (float2)(0.1f, 0.f), fseed);
        }
      }
    }

  if (val < 1.f)
  {
    val = val < 1.f ? sqrt(1.f - val) : 0.f;
    val *= (1.f - amplitude_random_ratio * amp_rnd);
  }
  else
  {
    val = 0.f;
  }

  return val;
}

void kernel hemisphere_field(global float *output,
                             global float *noise_x,
                             global float *noise_y,
                             global float *noise_distance,
                             global float *density_multiplier,
                             global float *size_multiplier,
                             const int     nx,
                             const int     ny,
                             float         kx,
                             float         ky,
                             const uint    seed,
                             const float   rmin,
                             const float   rmax,
                             const float   amplitude_random_ratio,
                             const float   density,
                             const float2  jitter,
                             const float   shift,
                             const int     has_noise_x,
                             const int     has_noise_y,
                             const int     has_noise_distance,
                             const int     has_density_multiplier,
                             const int     has_size_multiplier,
                             const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;
  float dr = has_noise_distance > 0 ? noise_distance[index] : 0.f;
  float dm = has_density_multiplier > 0 ? density_multiplier[index] : 1.f;
  float sm = has_size_multiplier > 0 ? size_multiplier[index] : 1.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  float val = hf_base(pos,
                      sm * rmin,
                      sm * rmax,
                      amplitude_random_ratio,
                      dm * density,
                      jitter,
                      shift + dr,
                      fseed);

  output[index] = val;
}

void kernel hemisphere_field_fbm(global float *output,
                                 global float *noise_x,
                                 global float *noise_y,
                                 global float *noise_distance,
                                 global float *density_multiplier,
                                 global float *size_multiplier,
                                 const int     nx,
                                 const int     ny,
                                 float         kx,
                                 float         ky,
                                 const uint    seed,
                                 const float   rmin,
                                 const float   rmax,
                                 const float   amplitude_random_ratio,
                                 const float   density,
                                 const float2  jitter,
                                 const float   shift,
                                 const int     octaves,
                                 const float   persistence,
                                 const float   lacunarity,
                                 const int     has_noise_x,
                                 const int     has_noise_y,
                                 const int     has_noise_distance,
                                 const int     has_density_multiplier,
                                 const int     has_size_multiplier,
                                 const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;
  float dr = has_noise_distance > 0 ? noise_distance[index] : 0.f;
  float dm = has_density_multiplier > 0 ? density_multiplier[index] : 1.f;
  float sm = has_size_multiplier > 0 ? size_multiplier[index] : 1.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  uint  seed_m = seed;
  for (int i = 0; i < octaves; i++)
  {
    uint  rng_state = wang_hash(seed_m++);
    float fseed = rand(&rng_state);

    float v = hf_base(nf * pos,
                      sm * rmin,
                      sm * rmax,
                      amplitude_random_ratio,
                      dm * density,
                      jitter,
                      shift + dr,
                      fseed);
    n = smax(n, v * na, 0.01f);
    na *= persistence;
    nf *= lacunarity;
  }

  output[index] = n;
}
)""
