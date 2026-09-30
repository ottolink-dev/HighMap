R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// --- BASE NOISE

// Wrap a lattice coordinate modulo a period so the noise tiles seamlessly.
// A non-positive period component disables wrapping on that axis (default,
// non-periodic behaviour). Lattice coordinates are integer-valued floats, so
// the result is exact at the period boundary.
float2 wrap_lattice(float2 c, const int2 period)
{
  if (period.x > 0) c.x = c.x - (float)period.x * floor(c.x / (float)period.x);
  if (period.y > 0) c.y = c.y - (float)period.y * floor(c.y / (float)period.y);
  return c;
}

// Scale a lattice period for the next fbm octave. Tiling stays exact only for
// integer lacunarity (the default is 2); the round keeps the period integral.
int2 scale_period(int2 period, const float lacunarity)
{
  if (period.x > 0) period.x = (int)(period.x * lacunarity + 0.5f);
  if (period.y > 0) period.y = (int)(period.y * lacunarity + 0.5f);
  return period;
}

float base_perlin(const float2 p, const float fseed, const int2 period)
{
  // lattice points
  float2 i = floor(p);
  float2 pi;
  float2 f = fract(p, &pi);

  // gradients at lattice corners (lattice indices wrapped for tileability)
  float2 g00 = grad22f(wrap_lattice(i + (float2)(0.0f, 0.0f), period), fseed);
  float2 g10 = grad22f(wrap_lattice(i + (float2)(1.0f, 0.0f), period), fseed);
  float2 g01 = grad22f(wrap_lattice(i + (float2)(0.0f, 1.0f), period), fseed);
  float2 g11 = grad22f(wrap_lattice(i + (float2)(1.0f, 1.0f), period), fseed);

  // offset vectors to corners
  float2 d00 = f - (float2)(0.0f, 0.0f);
  float2 d10 = f - (float2)(1.0f, 0.0f);
  float2 d01 = f - (float2)(0.0f, 1.0f);
  float2 d11 = f - (float2)(1.0f, 1.0f);

  // dot products between gradients and offset vectors
  float n00 = dot(g00, d00);
  float n10 = dot(g10, d10);
  float n01 = dot(g01, d01);
  float n11 = dot(g11, d11);

  // fade curve for interpolation
  float2 u = smoothstep5_f2(f);

  // bilinear interpolation
  float nx0 = lerp(n00, n10, u.x);
  float nx1 = lerp(n01, n11, u.x);
  float nxy = lerp(nx0, nx1, u.y);

  // return noise value
  return 1.42857f * nxy;
}

float base_simplex2(const float2 p, const float fseed)
{
  const float K1 = 0.366025403f; // (sqrt(3)-1)/2
  const float K2 = 0.211324865f; // (3-sqrt(3))/6

  // adjust to have consistent wavenb definition with other noises
  float2 p2 = 0.5f * p;

  // skew the input space to find the simplex cell
  float2 i = floor(p2 + dot(p2, (float2)(K1, K1)));
  float2 a = p2 - i + dot(i, (float2)(K2, K2));

  // determine which simplex we are in
  float  m = step(a.y, a.x);
  float2 o = (float2)(m, 1.0f - m);
  float2 b = a - o + (float2)(K2, K2);
  float2 c = a - 1.0f + 2.0f * (float2)(K2, K2);

  // compute gradients at each vertex
  float2 g0 = grad22f(i, fseed);
  float2 g1 = grad22f(i + o, fseed);
  float2 g2 = grad22f(i + (float2)(1.0f, 1.0f), fseed);

  // compute contributions from each vertex
  float n0 = max(0.5f - dot(a, a), 0.0f);
  float n1 = max(0.5f - dot(b, b), 0.0f);
  float n2 = max(0.5f - dot(c, c), 0.0f);

  n0 = n0 * n0 * n0 * n0 * dot(g0, a);
  n1 = n1 * n1 * n1 * n1 * dot(g1, b);
  n2 = n2 * n2 * n2 * n2 * dot(g2, c);

  // sum contributions and scale result (100.0 = 70.0 / 0.7)
  return 100.f * (n0 + n1 + n2);
}

float base_value(const float2 p, const float fseed, const int2 period)
{
  float2 i = floor(p);
  float2 pi;
  float2 f = fract(p, &pi);

  float v00 = hash12f(wrap_lattice(i + (float2)(0.0f, 0.0f), period), fseed);
  float v10 = hash12f(wrap_lattice(i + (float2)(1.0f, 0.0f), period), fseed);
  float v01 = hash12f(wrap_lattice(i + (float2)(0.0f, 1.0f), period), fseed);
  float v11 = hash12f(wrap_lattice(i + (float2)(1.0f, 1.0f), period), fseed);

  float2 u = smoothstep3_f2(f);

  float nx0 = lerp(v00, v10, u.x);
  float nx1 = lerp(v01, v11, u.x);
  float nxy = lerp(nx0, nx1, u.y);

  return 2.f * nxy - 1.f;
}

float base_value_cubic(const float2 p, const float fseed, const int2 period)
{
  float2 i = floor(p);
  float2 pi;
  float2 f = fract(p, &pi);

  // collect 16 lattice values around the point
  float v[4][4];
  for (int dy = -1; dy <= 2; dy++)
    for (int dx = -1; dx <= 2; dx++)
      v[dy + 1][dx + 1] = hash12f(wrap_lattice(i + (float2)(dx, dy), period),
                                  fseed);

  // cubic interpolation in the x-direction for each row
  float interp_row[4];
  for (int dy = 0; dy < 4; dy++)
    interp_row[dy] = cubic_interp(v[dy][0], v[dy][1], v[dy][2], v[dy][3], f.x);

  // cubic interpolation in the y-direction
  float value = cubic_interp(interp_row[0],
                             interp_row[1],
                             interp_row[2],
                             interp_row[3],
                             f.y);

  return 1.43f * (value - 0.5f);
}

float base_value_linear(const float2 p, const float fseed, const int2 period)
{
  float2 i = floor(p);
  float2 pi;
  float2 f = fract(p, &pi);

  float v00 = hash12f(wrap_lattice(i + (float2)(0.0f, 0.0f), period), fseed);
  float v10 = hash12f(wrap_lattice(i + (float2)(1.0f, 0.0f), period), fseed);
  float v01 = hash12f(wrap_lattice(i + (float2)(0.0f, 1.0f), period), fseed);
  float v11 = hash12f(wrap_lattice(i + (float2)(1.0f, 1.0f), period), fseed);

  float nx0 = lerp(v00, v10, f.x);
  float nx1 = lerp(v01, v11, f.x);
  float nxy = lerp(nx0, nx1, f.y);

  return 2.f * nxy - 1.f;
}

float base_worley(const float2 p, const float fseed, const int2 period)
{
  float2 i = floor(p);
  float2 pi;
  float2 f = fract(p, &pi);

  float min_dist = FLT_MAX;

  for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    {
      float2 dr = (float2)(dx, dy);
      float2 feature_point = dr + hash22f(wrap_lattice(i + dr, period), fseed);
      float2 diff = f - feature_point;
      float  dist = dot(diff, diff);

      min_dist = min(min_dist, dist);
    }

  // NB - squared distance
  return 1.66f * min_dist - 1.f;
}

// --- FBM

float base_perlin_fbm(const float2 p,
                      const int    octaves,
                      const float  weight,
                      const float  persistence,
                      const float  lacunarity,
                      const float  fseed,
                      const int2   period)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  int2  per = period;
  for (int i = 0; i < octaves; i++)
  {
    float v = base_perlin(p * nf, fseed, per);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }
  return n;
}

float base_perlin_billow_fbm(const float2 p,
                             const int    octaves,
                             const float  weight,
                             const float  persistence,
                             const float  lacunarity,
                             const float  fseed,
                             const int2   period)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  int2  per = period;
  for (int i = 0; i < octaves; i++)
  {
    float v = 2.f * fabs(base_perlin(p * nf, fseed, per)) - 1.f;
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }
  return n;
}

float base_perlin_half_fbm(const float2 p,
                           const int    octaves,
                           const float  weight,
                           const float  persistence,
                           const float  lacunarity,
                           const float  fseed,
                           const int2   period)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  int2  per = period;
  for (int i = 0; i < octaves; i++)
  {
    float v = max_smooth(base_perlin(p * nf, fseed, per), 0.f, 0.5f);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }
  return n;
}

float base_simplex2_fbm(const float2 p,
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
    float v = base_simplex2(p * nf, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_value_fbm(const float2 p,
                     const int    octaves,
                     const float  weight,
                     const float  persistence,
                     const float  lacunarity,
                     const float  fseed,
                     const int2   period)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  int2  per = period;
  for (int i = 0; i < octaves; i++)
  {
    float v = base_value(p * nf, fseed, per);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }
  return n;
}

float base_value_cubic_fbm(const float2 p,
                           const int    octaves,
                           const float  weight,
                           const float  persistence,
                           const float  lacunarity,
                           const float  fseed,
                           const int2   period)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  int2  per = period;
  for (int i = 0; i < octaves; i++)
  {
    float v = base_value_cubic(p * nf, fseed, per);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }
  return n;
}

float base_value_linear_fbm(const float2 p,
                            const int    octaves,
                            const float  weight,
                            const float  persistence,
                            const float  lacunarity,
                            const float  fseed,
                            const int2   period)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  int2  per = period;
  for (int i = 0; i < octaves; i++)
  {
    float v = base_value_linear(p * nf, fseed, per);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }
  return n;
}

float base_worley_fbm(const float2 p,
                      const int    octaves,
                      const float  weight,
                      const float  persistence,
                      const float  lacunarity,
                      const float  fseed,
                      const int2   period)
{
  float n = 0.f;
  float nf = 1.f;
  float na = 0.6f;
  int2  per = period;
  for (int i = 0; i < octaves; i++)
  {
    float v = base_worley(p * nf, fseed, per);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }
  return n;
}
)""
