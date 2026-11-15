R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

uint hash32(uint x)
{
  x ^= x >> 16;
  x *= 0x7feb352dU;
  x ^= x >> 15;
  x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

uint hash21u(uint x, uint y)
{
  uint h = x * 0x8da6b343u ^ y * 0xd8163841u;
  h ^= h >> 13;
  h *= 0x85ebca6bu;
  h ^= h >> 16;
  return h;
}

uint hash21u_seed(int ix, int iy, uint seed)
{
  uint h = (uint)ix * 0x8da6b343u ^ (uint)iy * 0xd8163841u ^ seed * 0x9e3779b9u;
  return hash32(h);
}

float hash12u_f(int ix, int iy, uint seed)
{
  return (float)(hash21u_seed(ix, iy, seed) >> 8) * (1.f / 16777216.f);
}

float2 hash22u_f(int ix, int iy, uint seed)
{
  uint h0 = hash21u_seed(ix, iy, seed);
  uint h1 = hash32(h0 ^ 0x9e3779b9u);
  return (float2)((float)(h0 >> 8) * (1.f / 16777216.f),
                  (float)(h1 >> 8) * (1.f / 16777216.f));
}

float2 grad22u_f(int ix, int iy, uint seed)
{
  float angle = 6.2831853f * hash12u_f(ix, iy, seed); // 2 * PI
  return (float2)(cos(angle), sin(angle));
}

float hash11f(float x, float fseed)
{
  uint h = hash32(as_uint(x) ^ as_uint(fseed) * 0x9e3779b9u);
  return (float)(h >> 8) * (1.f / 16777216.f);
}

float hash12f(float2 p, float fseed)
{
  uint h = hash32(as_uint(p.x) * 0x8da6b343u ^ as_uint(p.y) * 0xd8163841u ^
                  as_uint(fseed) * 0x9e3779b9u);
  return (float)(h >> 8) * (1.f / 16777216.f);
}

float hash13f(float3 p, float fseed)
{
  uint h = hash32(as_uint(p.x) * 0x8da6b343u ^ as_uint(p.y) * 0xd8163841u ^
                  as_uint(p.z) * 0x165667b1u ^ as_uint(fseed) * 0x9e3779b9u);
  return (float)(h >> 8) * (1.f / 16777216.f);
}

float2 hash22f(float2 p, float fseed)
{
  uint h0 = hash32(as_uint(p.x) * 0x8da6b343u ^ as_uint(p.y) * 0xd8163841u ^
                   as_uint(fseed) * 0x9e3779b9u);
  uint h1 = hash32(h0 ^ 0x9e3779b9u);
  return (float2)((float)(h0 >> 8) * (1.f / 16777216.f),
                  (float)(h1 >> 8) * (1.f / 16777216.f));
}

float2 hash22f_poly(float2 x, float fseed)
{
  const float2 k = {0.3183099f, 0.3678794f};
  x = x * k + k.yx;
  float2 qi;
  return fract(16.f * k * fract(x.x * x.y * (x.x + x.y) + fseed, &qi), &qi);
}

float3 hash33f(float3 p, float fseed)
{
  uint h0 = hash32(as_uint(p.x) * 0x8da6b343u ^ as_uint(p.y) * 0xd8163841u ^
                   as_uint(p.z) * 0x165667b1u ^ as_uint(fseed) * 0x9e3779b9u);
  uint h1 = hash32(h0 ^ 0x9e3779b9u);
  uint h2 = hash32(h1 ^ 0x85ebca6bu);
  return (float3)((float)(h0 >> 8) * (1.f / 16777216.f),
                  (float)(h1 >> 8) * (1.f / 16777216.f),
                  (float)(h2 >> 8) * (1.f / 16777216.f));
}

float2 grad22f(float2 p, float fseed)
{
  float angle = 6.2831853f * hash12f(p, fseed); // 2 * PI
  return (float2)(cos(angle), sin(angle));
}

uint pcg_hash(uint v)
{
  v = v * 747796405u + 2891336453u;
  uint word = ((v >> ((v >> 28u) + 4u)) ^ v) * 277803737u;
  return (word >> 22u) ^ word;
}

uint wang_hash(uint seed)
{
  seed = (seed ^ 61) ^ (seed >> 16);
  seed *= 9;
  seed = seed ^ (seed >> 4);
  seed *= 0x27d4eb2d;
  seed = seed ^ (seed >> 15);
  return seed;
}

void rand_xorshift(uint *p_rng_state)
{
  // Xorshift algorithm from George Marsaglia's paper
  *p_rng_state ^= (*p_rng_state << 13);
  *p_rng_state ^= (*p_rng_state >> 17);
  *p_rng_state ^= (*p_rng_state << 5);
}

float rand(uint *p_rng_state)
{
  rand_xorshift(p_rng_state);
  return *p_rng_state * (1.f / 4294967296.f);
}

void shuffle_array_i(int *v, int size, int seed)
{
  // use a simple linear congruential generator (LCG) for
  // pseudo-random numbers
  int rng = seed;

  // fisher-Yates shuffle
  for (int i = size - 1; i > 0; i--)
  {
    // generate a random index between 0 and i
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF; // LCG
    int j = rng % (i + 1);

    int temp_v = v[i];
    v[i] = v[j];
    v[j] = temp_v;
  }
}

void shuffle_array_f(float *v, int size, int seed)
{
  // use a simple linear congruential generator (LCG) for
  // pseudo-random numbers
  int rng = seed;

  // fisher-Yates shuffle
  for (int i = size - 1; i > 0; i--)
  {
    // generate a random index between 0 and i
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF; // LCG
    int j = rng % (i + 1);

    float temp_v = v[i];
    v[i] = v[j];
    v[j] = temp_v;
  }
}
)""
