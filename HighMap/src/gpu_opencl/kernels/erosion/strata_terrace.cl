R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
float stratify_rand(const float val,
                    const float shift,
                    const float kz,
                    const float gamma,
                    const bool  linear_gamma,
                    const float gamma_noise_ratio,
                    const float fseed)
{
  float out = val + shift;
  out *= kz;

  // apply recruve function using a cell-based (but pointwise)
  // approach
  float i = floor(out);
  float u = out - i;

  float rnd = 2.f * hash11f(i, fseed) - 1.f;
  float gamma_r = gamma * (1.f + gamma_noise_ratio * rnd);
  gamma_r = clamp(gamma_r, 0.01f, 10.f);

  if (linear_gamma)
  {
    // sharp linear gamma
    float a = pow(1.f / gamma_r, 1.f / (gamma_r - 1.f));
    float b = pow(gamma_r, -gamma_r / (gamma_r - 1.f));

    if (u < a)
      u *= b / a;
    else
      u = b + (1.f - b) * (u - a) / (1.f - a);
  }
  else
  {
    // smooth gamma
    float ce = 50.f / gamma_r;
    u = pow(u, gamma_r) * (1.f - exp(-ce * u));
  }
  // revert back value to input range
  out = (i + u) / kz - shift;

  return out;
}

void kernel strata_terrace(global float *output,
                           global float *noise,
                           const int     nx,
                           const int     ny,
                           const float   gamma,
                           const uint    seed,
                           const int     use_linear_gamma,
                           const float   kz,
                           const float   gamma_noise_ratio,
                           const float   slope,
                           const float   angle,
                           const int     has_noise,
                           const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int    index = linear_index(g.x, g.y, nx);
  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, 0.f, 0.f, bbox);

  // parameters
  const bool linear_gamma = use_linear_gamma == 0 ? false : true;
  float      shift = has_noise > 0 ? noise[index] : 0.f;

  float  alpha = angle / 180.f * 3.141592f;
  float2 dir = (float2)(cos(alpha), sin(alpha));
  float  dr = dot(pos, dir);
  shift += slope * dr;

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  // correct wavenumber with slope to avoid very high-frequencies for
  // high slopes
  float kz_corrected = kz / (1.f + slope);

  float val = output[index];
  val = stratify_rand(val,
                      shift,
                      kz_corrected,
                      gamma,
                      linear_gamma,
                      gamma_noise_ratio,
                      fseed);

  // out
  output[index] = val;
}
)""
