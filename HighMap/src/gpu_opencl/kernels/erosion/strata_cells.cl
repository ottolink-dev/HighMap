R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
float base_strata_cells(const float  val,
                        const float2 p,
                        const float  gamma,
                        const float  gamma_lateral,
                        const float  amp,
                        const float  noise_amp,
                        const bool   absolute_displacement,
                        const float  occurence_probability,
                        const float  fseed)
{
  float2 i = floor(p);
  float2 pi;
  float2 f = fract(p, &pi);

  float out = val;
  float rnd01 = hash12f(i, fseed);
  float alpha = rnd01 * 90.f / 180.f * 3.1416f;

  if (alpha != 0.f)
  {
    float ca = cos(alpha);
    float sa = sin(alpha);
    f -= 0.5f;
    f = (float2)(f.x * ca - f.y * sa, f.x * sa + f.y * ca);
    f += 0.5f;

    f.x = smin(f.x, 1.f, 0.1f);
    f.x = smax(f.x, 0.f, 0.1f);
    f.y = smin(f.y, 1.f, 0.1f);
    f.y = smax(f.y, 0.f, 0.1f);
  }

  // skip some cells arbitrarily
  if (rnd01 > occurence_probability) return val;

  float rnd = 2.f * rnd01 - 1.f;

  // transverse profile
  float ay = (4.f * f.y * (1.f - f.y));
  ay = pow(ay, gamma_lateral);

  // axial profile
  // (1.f + noise_amp * rnd)

  float a = pow(gamma, gamma) / pow(gamma + 1.f, gamma + 1.f);
  float u = pow(f.x, gamma) * (1.f - f.x) / a;
  u = clamp(u - amp * rnd01, 0.f, 1.f);

  if (!absolute_displacement) u *= val;

  out = val + u * amp * (1.f + noise_amp * rnd) * ay;

  return out;
}

void kernel strata_cells(global float *output,
                         global float *mask,
                         global float *noise_x,
                         global float *noise_y,
                         const int     nx,
                         const int     ny,
                         const float   kx,
                         const float   ky,
                         const uint    seed,
                         const float   gamma,
                         const float   gamma_lateral,
                         const float   angle,
                         const float   amp,
                         const float   noise_amp,
                         const int     absolute_displacement,
                         const float   occurence_probability,
                         const int     has_mask,
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

  // add rotation

  float alpha = angle / 180.f * 3.14159f;
  float ca = cos(alpha);
  float sa = sin(alpha);
  pos = (float2)(pos.x * ca - pos.y * sa, pos.x * sa + pos.y * ca);

  const float val0 = output[index];
  float       val = base_strata_cells(val0,
                                pos,
                                gamma,
                                gamma_lateral,
                                amp,
                                noise_amp,
                                absolute_displacement,
                                occurence_probability,
                                fseed);

  // input mask
  if (has_mask)
  {
    float t = mask[index];
    val = lerp(val0, val, t);
  }

  // out
  output[index] = val;
}
)""
