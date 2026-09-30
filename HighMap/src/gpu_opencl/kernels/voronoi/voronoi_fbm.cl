R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

float base_voronoi_f1_fbm(const float2 p,
                          const float2 jitter,
                          const float  k_smoothing,
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
    float v = base_voronoi_f1(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_f2_fbm(const float2 p,
                          const float2 jitter,
                          const float  k_smoothing,
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
    float v = base_voronoi_f2(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_f1tf2_fbm(const float2 p,
                             const float2 jitter,
                             const float  k_smoothing,
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
    float v = base_voronoi_f1tf2(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_f1df2_fbm(const float2 p,
                             const float2 jitter,
                             const float  k_smoothing,
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
    float v = base_voronoi_f1df2(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_f2mf1_fbm(const float2 p,
                             const float2 jitter,
                             const float  k_smoothing,
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
    float v = base_voronoi_f2mf1(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_edge_distance_fbm(const float2 p,
                                     const float2 jitter,
                                     const float  k_smoothing,
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
    float v = base_voronoi_edge_distance(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_edge_distance_exp_fbm(const float2 p,
                                         const float2 jitter,
                                         const float  k_smoothing,
                                         const float  exp_sigma,
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
    float v = base_voronoi_edge_distance_exp(p * nf,
                                             jitter,
                                             k_smoothing,
                                             exp_sigma,
                                             fseed);
    n = max(n, v * na);
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_constant_fbm(const float2 p,
                                const float2 jitter,
                                const float  k_smoothing,
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
    float v = base_voronoi_constant(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}

float base_voronoi_constant_f2mf1_fbm(const float2 p,
                                      const float2 jitter,
                                      const float  k_smoothing,
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
    float v = base_voronoi_constant_f2mf1(p * nf, jitter, k_smoothing, fseed);
    n += v * na;
    na *= (1.f - weight) + weight * min(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }
  return n;
}
)""
