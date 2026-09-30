R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

/* Clean-room port of the flow-coupled particle erosion model from
 * erosiv/soillib (LGPL-3, Nicholas McDonald): persistent discharge/momentum
 * fields couple particles through the mean flow; a bank-stability debris
 * flow runs in the same solver loop against a separate sediment layer. */

/* one deterministic uniform in [0,1) per (seed, step, sample, salt) */
float mcd_uniform(const uint seed,
                  const uint step,
                  const uint n,
                  const uint salt)
{
  return (float)(hash21u(n ^ (step * 0x9e3779b9u), seed + salt) >> 8) /
         16777216.f;
}

/* 5-point one-sided-aware gradient of (bed + sed), port of soillib
 * gather.hpp lerp5_t::gather + grad. Returns d(height)/d(cell) scaled by
 * z_m / cell_m (grad's s = scale.z / scale.xy). */
float2 mcd_grad5(global const float *bed,
                 global const float *sed,
                 const int2          p,
                 const int           nx,
                 const int           ny,
                 const float         zs) /* zs = z_m / cell_m */
{
  float vx[5];
  int   ox[5];
  float vy[5];
  int   oy[5];

  for (int k = 0; k < 5; ++k)
  {
    int xi = p.x - 2 + k;
    int yj = p.y - 2 + k;

    ox[k] = (xi < 0 || xi >= nx || p.y < 0 || p.y >= ny) ? 1 : 0;
    if (!ox[k])
    {
      int idx = linear_index(xi, p.y, nx);
      vx[k] = bed[idx] + sed[idx];
    }
    else
      vx[k] = 0.f;

    oy[k] = (yj < 0 || yj >= ny || p.x < 0 || p.x >= nx) ? 1 : 0;
    if (!oy[k])
    {
      int idx = linear_index(p.x, yj, nx);
      vy[k] = bed[idx] + sed[idx];
    }
    else
      vy[k] = 0.f;
  }

  float gx = 0.f;
  if (!ox[0] && !ox[4])
    gx = (1.f * vx[0] - 8.f * vx[1] + 8.f * vx[3] - 1.f * vx[4]) / 12.f;
  else if (!ox[0] && !ox[3])
    gx = (1.f * vx[0] - 6.f * vx[1] + 3.f * vx[2] + 2.f * vx[3]) / 6.f;
  else if (!ox[0] && !ox[2])
    gx = (1.f * vx[0] - 4.f * vx[1] + 3.f * vx[2]) / 2.f;
  else if (!ox[1] && !ox[4])
    gx = (-2.f * vx[1] - 3.f * vx[2] + 6.f * vx[3] - 1.f * vx[4]) / 6.f;
  else if (!ox[2] && !ox[4])
    gx = (-3.f * vx[2] + 4.f * vx[3] - 1.f * vx[4]) / 2.f;
  else if (!ox[1] && !ox[3])
    gx = (-1.f * vx[1] + 1.f * vx[3]) / 2.f;
  else if (!ox[2] && !ox[3])
    gx = (-1.f * vx[2] + 1.f * vx[3]) / 1.f;
  else if (!ox[1] && !ox[2])
    gx = (-1.f * vx[1] + 1.f * vx[2]) / 1.f;

  float gy = 0.f;
  if (!oy[0] && !oy[4])
    gy = (1.f * vy[0] - 8.f * vy[1] + 8.f * vy[3] - 1.f * vy[4]) / 12.f;
  else if (!oy[0] && !oy[3])
    gy = (1.f * vy[0] - 6.f * vy[1] + 3.f * vy[2] + 2.f * vy[3]) / 6.f;
  else if (!oy[0] && !oy[2])
    gy = (1.f * vy[0] - 4.f * vy[1] + 3.f * vy[2]) / 2.f;
  else if (!oy[1] && !oy[4])
    gy = (-2.f * vy[1] - 3.f * vy[2] + 6.f * vy[3] - 1.f * vy[4]) / 6.f;
  else if (!oy[2] && !oy[4])
    gy = (-3.f * vy[2] + 4.f * vy[3] - 1.f * vy[4]) / 2.f;
  else if (!oy[1] && !oy[3])
    gy = (-1.f * vy[1] + 1.f * vy[3]) / 2.f;
  else if (!oy[2] && !oy[3])
    gy = (-1.f * vy[2] + 1.f * vy[3]) / 1.f;
  else if (!oy[1] && !oy[2])
    gy = (-1.f * vy[1] + 1.f * vy[2]) / 1.f;

  return (float2)(gx * zs, gy * zs);
}

void kernel mcdonald_solve(global float       *bed,
                           global float       *sed,
                           global const float *dis,
                           global const float *mx,
                           global const float *my,
                           global float       *tr_d,
                           global float       *tr_mx,
                           global float       *tr_my,
                           global float       *moisture_map,
                           const int           nx,
                           const int           ny,
                           const int           n_samples,
                           const uint          seed,
                           const uint          step,
                           const float         cell_m,
                           const float         z_m,
                           const float         time_step,
                           const float         rainfall,
                           const float         evap_rate,
                           const float         gravity,
                           const float         viscosity,
                           const float         bed_shear,
                           const float         deposition_rate,
                           const float         suspension_rate,
                           const float         exit_slope,
                           const int           maxage,
                           const int           has_moisture_map)
{
  int n = get_global_id(0);
  if (n >= n_samples) return;

  const float2 cl = (float2)(cell_m, cell_m);
  const float  cl_len = length(cl);
  const float  Ac = cell_m * cell_m;
  const float  zs = z_m / cell_m;

  float2 pos = (float2)(mcd_uniform(seed, step, (uint)n, 0u) * (float)nx,
                        mcd_uniform(seed, step, (uint)n, 1u) * (float)ny);
  int    ci = (int)pos.x;
  int    cj = (int)pos.y;
  if (ci >= nx || cj >= ny) return;
  int find = linear_index(ci, cj, nx);

  const float P = 1.f / ((float)nx * (float)ny);
  const float Z = Ac * z_m;
  const float Q = P * (float)n_samples;

  float m = (has_moisture_map != 0) ? moisture_map[find] : 1.f;
  float vol = Ac * rainfall * m;
  float sedm = 0.f;

  float2 grad = mcd_grad5(bed, sed, (int2)(ci, cj), nx, ny, zs);
  float3 normal = normalize((float3)(-grad.x, -grad.y, 1.f));

  float2 avg = (float2)(0.f, 0.f);
  if (dis[find] > 0.f) avg = (float2)(mx[find], my[find]) / dis[find];

  float2 speed = gravity * (float2)(normal.x, normal.y) + viscosity * avg;
  if (length(speed) == 0.f) return;

  float  ds = cl_len / length(speed);
  float2 npos = pos + ds * (speed / cl);
  float2 dspeed = speed;

  for (int age = 0; age < maxage; ++age)
  {
    atomic_add_float(&tr_d[find], (1.f / P / (float)n_samples) * vol);
    atomic_add_float(&tr_mx[find],
                     (1.f / P / (float)n_samples) * vol * dspeed.x);
    atomic_add_float(&tr_my[find],
                     (1.f / P / (float)n_samples) * vol * dspeed.y);

    float discharge = dis[find];
    float slope = -exit_slope;
    float h0 = (bed[find] + sed[find]) * z_m;
    float h1 = h0 + slope * cl_len;

    int ni = (int)npos.x;
    int nj = (int)npos.y;
    if (ni >= 0 && ni < nx && nj >= 0 && nj < ny)
    {
      int nind = linear_index(ni, nj, nx);
      h1 = (bed[nind] + sed[nind]) * z_m;
      slope = (h1 - h0) / cl_len;
    }

    float alpha = (slope < 0.f) ? 1.f : 0.f;
    float suspend = time_step * suspension_rate * vol * slope * alpha *
                    pow(discharge, 0.4f);
    float deposit = time_step * deposition_rate * sedm;

    /* multi-material mass transfer (reference "Multi-Material" branch);
     * clamp written as fmin(fmax(...)) to match glm componentwise behavior */
    float transfer = deposit + suspend;
    float maxtransfer = 0.1f * slope * cl_len / z_m * Z * Q;
    float tmin = transfer * fmin(1.f, fabs(maxtransfer / suspend));
    float tmax = sedm;
    transfer = fmin(fmax(transfer, tmin), tmax);

    // clamp erosion transfer to prevent numerical explosion
    float max_erode = 0.5f * fmax(0.f, h0) * Ac * Q;
    transfer = fmax(transfer, -max_erode);

    if (transfer > 0.f)
    {
      atomic_add_float(&sed[find], transfer / Z / Q);
      sedm -= transfer;
    }
    else if (transfer < 0.f)
    {
      float maxtransfer2 = 0.1f * sed[find] * Z * Q;
      float t1 = transfer * fmin(1.f, fabs(maxtransfer2 / transfer));
      atomic_add_float(&sed[find], t1 / Z / Q);
      sedm -= t1;

      transfer -= t1;
      atomic_add_float(&bed[find], transfer / Z / Q);
      sedm -= transfer;
    }

    vol = 1.f / (1.f + ds * evap_rate) * vol;

    pos = npos;
    ci = (int)pos.x;
    cj = (int)pos.y;
    if (ci < 0 || ci >= nx || cj < 0 || cj >= ny) break;
    find = linear_index(ci, cj, nx);

    grad = mcd_grad5(bed, sed, (int2)(ci, cj), nx, ny, zs);
    normal = normalize((float3)(-grad.x, -grad.y, 1.f));

    float2 avg2 = (float2)(0.f, 0.f);
    if (dis[find] > 0.f) avg2 = (float2)(mx[find], my[find]) / dis[find];

    speed = speed + ds * gravity * (float2)(normal.x, normal.y);

    float k12 = bed_shear + viscosity;
    speed = 1.f / (1.f + ds * k12) * speed +
            ds * viscosity / (1.f + ds * k12) * avg2;
    dspeed = 1.f / (1.f + ds * k12) * dspeed;

    if (length(speed) == 0.f) break;

    ds = cl_len / length(speed);
    npos = pos + ds * (speed / cl);
  }
}

void kernel mcdonald_filter(global float       *dis,
                            global float       *mx,
                            global float       *my,
                            global const float *tr_d,
                            global const float *tr_mx,
                            global const float *tr_my,
                            const int           nx,
                            const int           ny,
                            const float         lrate)
{
  int2 g = {get_global_id(0), get_global_id(1)};
  if (g.x >= nx || g.y >= ny) return;
  int idx = linear_index(g.x, g.y, nx);

  dis[idx] = mix(dis[idx], tr_d[idx], lrate);
  mx[idx] = mix(mx[idx], tr_mx[idx], lrate);
  my[idx] = mix(my[idx], tr_my[idx], lrate);
}

void kernel mcdonald_debris(global float *bed,
                            global float *sed,
                            const int     nx,
                            const int     ny,
                            const int     n_samples,
                            const uint    seed,
                            const uint    step,
                            const float   cell_m,
                            const float   z_m,
                            const float   time_step,
                            const float   gravity,
                            const float   crit_slope,
                            const float   settle_rate,
                            const float   thermal_rate)
{
  int n = get_global_id(0);
  if (n >= n_samples)
    return; /* reference guards ind >= elem; identical while samples <= cells */

  const float2 cl = (float2)(cell_m, cell_m);
  const float  Ac = cell_m * cell_m;
  const float  zs = z_m / cell_m;
  const float  P = 1.f / ((float)nx * (float)ny);
  const float  Q = P * (float)n_samples;

  float mass = 0.f;

  float2 pos = (float2)(mcd_uniform(seed, step, (uint)n, 2u) * (float)nx,
                        mcd_uniform(seed, step, (uint)n, 3u) * (float)ny);

  for (int age = 0; age < 256; ++age)
  {
    int pi = (int)pos.x;
    int pj = (int)pos.y;
    if (pi < 0 || pi >= nx || pj < 0 || pj >= ny) return;

    float2 grad = mcd_grad5(bed, sed, (int2)(pi, pj), nx, ny, zs);
    float3 normal = normalize((float3)(-grad.x, -grad.y, 1.f));
    float2 speed = gravity * (float2)(normal.x, normal.y);

    float2 npos = pos;
    if (length(speed) > 0.f) npos = pos + normalize(speed);

    int qi = (int)npos.x;
    int qj = (int)npos.y;
    if (qi < 0 || qi >= nx || qj < 0 || qj >= ny) return;

    int find = linear_index(pi, pj, nx);
    int nind = linear_index(qi, qj, nx);

    float dist = length(cl * (npos - pos));
    pos = npos;

    float hf = z_m * bed[find] + fmax(0.f, z_m * sed[find]);
    float hf1 = fmax(0.f, z_m * sed[find]);
    float hn = z_m * bed[nind] + fmax(0.f, z_m * sed[nind]);

    float stable1 = hn + crit_slope * dist;

    float deposit = time_step * settle_rate * mass;
    float suspend = -time_step * thermal_rate * fmax(0.f, hf - stable1) * Ac;
    float transfer = deposit + suspend;
    if (transfer == 0.f) continue;

    if (transfer > 0.f)
    {
      float maxtransfer = fmax(0.f, stable1 - hf) * Ac * Q;
      transfer = fmin(transfer, maxtransfer);
      transfer = fmin(transfer, mass);
      transfer = fmax(0.f, transfer);

      atomic_add_float(&sed[find], transfer / Q / z_m / Ac);
      mass -= transfer;
    }
    else
    {
      float maxtransfer = fmax(0.f, hf - stable1) * Ac * Q;
      transfer = -fmin(-transfer, maxtransfer);

      // clamp erosion transfer to prevent numerical explosion
      float max_erode = 0.5f * fmax(0.f, hf) * Ac * Q;
      transfer = fmax(transfer, -max_erode);

      float maxt1 = hf1 * Ac * Q;
      float t1 = transfer * fmin(1.f, fabs(maxt1 / transfer));
      atomic_add_float(&sed[find], t1 / Q / z_m / Ac);
      mass -= t1;

      transfer -= t1;
      atomic_add_float(&bed[find], transfer / Q / Ac / z_m);
      mass -= transfer;
    }
  }
}
)""
