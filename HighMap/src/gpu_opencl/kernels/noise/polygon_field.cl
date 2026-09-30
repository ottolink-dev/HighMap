R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

inline float sdf_segment(float2 p, float2 a, float2 b)
{
  float2 pa = p - a;
  float2 ba = b - a;
  float  h = clamp(dot(pa, ba) / dot(ba, ba), 0.f, 1.f);
  float2 proj = a + ba * h;
  return length(p - proj);
}

inline float sdf_polygon(float2 p, const float2 *verts, int n)
{
  float d = FLT_MAX;
  int   inside = 0;

  for (int i = 0, j = n - 1; i < n; j = i++)
  {
    float2 a = verts[j];
    float2 b = verts[i];

    d = fmin(d, sdf_segment(p, a, b));

    int cond = ((a.y > p.y) != (b.y > p.y));
    if (cond)
    {
      float x = (b.x - a.x) * (p.y - a.y) / (b.y - a.y) + a.x;
      if (p.x < x) inside = !inside;
    }
  }

  return inside ? -d : d;
}

float pf_base(const float2 pos,
              const float  rmin,
              const float  rmax,
              const float  clamping_dist,
              const float  clamping_k,
              const int    n_vertices_min,
              const int    n_vertices_max,
              const float  density,
              const float2 jitter,
              const float  shift,
              const float  fseed)
{
  const int n = 64; // max polygon vertices
  float2    verts[64];

  int nvmax = min(n_vertices_max, n); // failsafe

  float  val = 0.f;
  float2 pos_i = floor(pos);

  for (int dx = -2; dx <= 2; dx++)
    for (int dy = -2; dy <= 2; dy++)
    {
      float2 pos_nbrs = pos_i + (float2)(dx, dy);

      // occurence probability
      float does_spawn_proba = hash12f(pos_nbrs, fseed);

      if (does_spawn_proba < density)
      {
        // nb of vertices
        float rnd_vert = hash12f(pos_nbrs, fseed);
        int   n_vertices = rint(n_vertices_min +
                              rnd_vert * (nvmax - n_vertices_min));

        // polar random generation
        for (int k = 0; k < n_vertices; k++)
        {
          float angle = (2.f * M_PI_F * k) / n_vertices;
          float radius = rmin + (rmax - rmin) *
                                    hash12f(pos_nbrs + (float2)(k, 0), fseed);
          float2 delta = jitter * hash22f(pos_nbrs, fseed);

          verts[k] = pos_nbrs + delta +
                     (float2)(cos(angle), sin(angle)) * radius;
        }

        // keep inner distance only
        val += max(0.f, -(sdf_polygon(pos, verts, n_vertices) - shift));
      }
    }

  return smin(val, clamping_dist, clamping_k);
}

void kernel polygon_field(global float *output,
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
                          const float   clamping_dist,
                          const float   clamping_k,
                          const int     n_vertices_min,
                          const int     n_vertices_max,
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

  float val = pf_base(pos,
                      sm * rmin,
                      sm * rmax,
                      clamping_dist,
                      clamping_k,
                      n_vertices_min,
                      n_vertices_max,
                      dm * density,
                      jitter,
                      shift + dr,
                      fseed);

  output[index] = val;
}

void kernel polygon_field_fbm(global float *output,
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
                              const float   clamping_dist,
                              const float   clamping_k,
                              const int     n_vertices_min,
                              const int     n_vertices_max,
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

    float v = pf_base(nf * pos,
                      sm * rmin,
                      sm * rmax,
                      clamping_dist,
                      clamping_k,
                      n_vertices_min,
                      n_vertices_max,
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
