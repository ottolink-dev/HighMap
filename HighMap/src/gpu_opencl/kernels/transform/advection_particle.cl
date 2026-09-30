R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel advection_particle(global float *field,
                               global float *dx,
                               global float *dy,
                               global float *out,
                               global float *count,
                               global float *advection_mask,
                               const int     nx,
                               const int     ny,
                               const int     nparticles,
                               const uint    seed,
                               const float   advection_sign,
                               const float   advection_length,
                               const float   value_persistence,
                               const float   inertia,
                               const int     has_advection_mask,
                               const int     num_fields)
{
  int id = get_global_id(0); // particle id
  if (id >= nparticles) return;

  uint rng_state = pcg_hash(seed + id * 2u);

  const float drx = 1.f;
  const float dry = 1.f;
  const float eps = 1e-9f;
  const float dist_max = advection_length * nx;
  float       dist = 0.f;

  float2 pos;
  if (has_advection_mask == 0)
  {
    pos = (float2)((nx - 1) * rand(&rng_state), (ny - 1) * rand(&rng_state));
  }
  else
  {
    // find a source point outside the advection zone
    int tries_count = 0;
    while (tries_count < 100)
    {
      tries_count++;
      pos = (float2)((nx - 1) * rand(&rng_state), (ny - 1) * rand(&rng_state));

      int2 g_try = (int2)((int)pos.x, (int)pos.y);
      if (advection_mask[linear_index_g(g_try, nx)] > 0) break;
    }
  }

  int2 g = (int2)((int)pos.x, (int)pos.y);
  int  idx = linear_index_g(g, nx);

  float val[64];
  int   num_f = num_fields;
  if (num_f > 64) num_f = 64;
  int stride = nx * ny;

  for (int f = 0; f < num_f; ++f)
  {
    val[f] = field[f * stride + idx];
  }

  float2 dir = (float2)(0.f, 0.f);
  float2 dir_prev = (float2)(0.f, 0.f);

  while (dist < dist_max)
  {
    // particle direction
    float2 vel = (float2)(dx[idx], dy[idx]);

    float len = length(vel);
    if (len > eps)
    {
      // (NB - if too flat, keep the direction of the previous step)
      dir = vel / len;
      dir = inertia * dir_prev + (1.f - inertia) * dir;
      dir = normalize(dir);
      dir_prev = dir;
    }

    dist += 1.f;

    // move
    float2 pos_prev = pos;

    pos.x += advection_sign * dir.x * drx;
    pos.y += advection_sign * dir.y * dry;

    if (fabs(pos.x - pos_prev.x) < eps && fabs(pos.y - pos_prev.y) < eps)
      continue;

    g = (int2)((int)pos.x, (int)pos.y);
    idx = linear_index_g(g, nx);

    // stop if the particle reaches the domain limits
    if (g.x < 0 || g.x > nx - 1 || g.y < 0 || g.y > ny - 1) break;

    // stop if outside the advection mask
    if (has_advection_mask > 0)
      if (advection_mask[idx] == 0) break;

    if (true)
    {
      for (int f = 0; f < num_f; ++f)
      {
        val[f] = lerp(field[f * stride + idx], val[f], value_persistence);
        out[f * stride + idx] += val[f];
      }
      count[idx] += 1.f;
    }
    else
    {
      int ir = 2;
      if (g.x < ir || g.x > nx - 1 - ir || g.y < ir || g.y > ny - 1 - ir) break;

      for (int f = 0; f < num_f; ++f)
      {
        val[f] = lerp(field[f * stride + idx], val[f], value_persistence);
      }

      float sum = 0.f;
      for (int p = -ir; p < ir; ++p)
        for (int q = -ir; q < ir; ++q)
        {
          float w = sqrt((float)(ir - abs(p)) * (ir - abs(q)));
          sum += w;
        }

      for (int p = -ir; p < ir; ++p)
        for (int q = -ir; q < ir; ++q)
        {
          float w = sqrt((float)(ir - abs(p)) * (ir - abs(q))) / sum;
          int2  gd = (int2)(g.x + p, g.y + q);
          int   out_idx = linear_index_g(gd, nx);
          for (int f = 0; f < num_f; ++f)
          {
            out[f * stride + out_idx] += val[f] * w;
          }
          count[out_idx] += w;
        }
    }
  }
}
)""
