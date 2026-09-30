R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// https://iquilezles.org/articles/voronoilines/
float voronoi_edge_distance_fct(const float2 x,
                                const float2 jitter,
                                const float  fseed)
{
  float2 p = floor(x);
  float2 pi;
  float2 f = fract(x, &pi);

  float2 mb;
  float2 mr;
  float2 df = (float2)(0.1f, 0.1f);

  float res = 8.f;
  for (int j = -1; j <= 1; j++)
    for (int i = -1; i <= 1; i++)
    {
      float2 b = (float2)(i, j);
      float2 r = b - f +
                 jitter * (float2)(hash12f(p + b, fseed),
                                   hash12f(p + b + df, fseed));
      float d = dot(r, r);

      if (d < res)
      {
        res = d;
        mr = r;
        mb = b;
      }
    }

  res = 8.f;
  for (int j = -2; j <= 2; j++)
    for (int i = -2; i <= 2; i++)
    {
      float2 b = mb + (float2)(i, j);
      float2 r = b - f +
                 jitter * (float2)(hash12f(p + b, fseed),
                                   hash12f(p + b + df, fseed));
      if (dot(mr - r, mr - r) > FLT_MIN)
      {
        float d = dot(0.5f * (mr + r), normalize(r - mr));
        res = min(res, d);
      }
    }

  return res;
}

void kernel voronoi_edge_distance(global float *output,
                                  global float *ctrl_param,
                                  global float *noise_x,
                                  global float *noise_y,
                                  const int     nx,
                                  const int     ny,
                                  const float   kx,
                                  const float   ky,
                                  const uint    seed,
                                  const float2  jitter,
                                  const int     has_ctrl_param,
                                  const int     has_noise_x,
                                  const int     has_noise_y,
                                  const float4  bbox)
{
  float width_ratio = 0.1f;

  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float ct = has_ctrl_param > 0 ? ctrl_param[index] : 1.f;
  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  float d = voronoi_edge_distance_fct(pos, ct * jitter, fseed);

  output[index] = d;
}
)""
