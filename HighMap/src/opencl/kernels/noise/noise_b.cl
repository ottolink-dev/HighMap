R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// --- OCL KERNEL

float2 apply_stretch(float2 p, float cos_a, float sin_a, float stretch_k)
{
  float inv_stretch = 1.f / fmax(stretch_k, 1e-6f);
  float xl = p.x * cos_a + p.y * sin_a;
  float yl = -p.x * sin_a + p.y * cos_a;
  xl *= inv_stretch;
  return (float2)(xl * cos_a - yl * sin_a, xl * sin_a + yl * cos_a);
}

float eval_base_noise(const int    noise_id,
                      const float2 pos,
                      const uint   seed,
                      const int2   period)
{
  if (noise_id == 1)
  {
    return base_perlin(pos, seed, period);
  }
  else if (noise_id == 2)
  {
    float value = base_perlin(pos, seed, period);
    return 2.f * fabs(value) - 1.f;
  }
  else if (noise_id == 3)
  {
    float value = base_perlin(pos, seed, period);
    return max_smooth(value, 0.f, 0.5f);
  }
  else if (noise_id == 4 || noise_id == 5)
  {
    return base_simplex2(pos, seed);
  }
  else if (noise_id == 6)
  {
    return base_value(pos, seed, period);
  }
  else if (noise_id == 7)
  {
    return base_value_cubic(pos, seed, period);
  }
  else if (noise_id == 9)
  {
    return base_value_linear(pos, seed, period);
  }
  else if (noise_id == 10)
  {
    return base_worley(pos, seed, period);
  }
  else if (noise_id == 11)
  {
    return base_worley_double(pos, seed, 0.5f, 0.f, period);
  }
  else if (noise_id == 12)
  {
    return base_worley_value(pos, seed, period);
  }
  return 0.f;
}

void kernel noise(global float *output,
                  global float *noise_x,
                  global float *noise_y,
                  const int     nx,
                  const int     ny,
                  const int     noise_id,
                  const float   kx,
                  const float   ky,
                  const uint    seed,
                  const int     has_noise_x,
                  const int     has_noise_y,
                  const int     period_x,
                  const int     period_y,
                  const float4  bbox,
                  const float   angle,
                  const float   stretch)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  int2 period = {period_x, period_y};

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  if (stretch != 1.f)
  {
    float theta = angle * 0.017453292519943295f;
    pos = apply_stretch(pos, cos(theta), sin(theta), stretch);
  }

  output[index] = eval_base_noise(noise_id, pos, seed, period);
}

void kernel noise_fbm(global float *output,
                      global float *ctrl_param,
                      global float *noise_x,
                      global float *noise_y,
                      const int     nx,
                      const int     ny,
                      const int     noise_id,
                      float         kx,
                      float         ky,
                      const uint    seed,
                      const int     octaves,
                      const float   weight,
                      const float   persistence,
                      const float   lacunarity,
                      const int     has_ctrl_param,
                      const int     has_noise_x,
                      const int     has_noise_y,
                      const int     period_x,
                      const int     period_y,
                      const float4  bbox,
                      const float   angle,
                      const float   stretch,
                      const float   stretch_decay)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float ct = has_ctrl_param > 0 ? ctrl_param[index] : 1.f;
  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  int2 period = {period_x, period_y};

  float  new_weight = (1.f - ct) + ct * weight;
  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  if (stretch != 1.f)
  {
    float theta = angle * 0.017453292519943295f;
    float cos_a = cos(theta);
    float sin_a = sin(theta);
    float decay_base = 1.f - clamp(stretch_decay, 0.f, 1.f);

    float n = 0.f;
    float nf = 1.f;
    float na = 0.6f;
    int2  per = period;

    for (int i = 0; i < octaves; i++)
    {
      float  stretch_k = 1.f + (stretch - 1.f) * pow(decay_base, (float)i);
      float2 p_s = apply_stretch(pos * nf, cos_a, sin_a, stretch_k);
      float  v = eval_base_noise(noise_id, p_s, seed, per);
      n += v * na;
      na *= (1.f - new_weight) + new_weight * min(v + 1.f, 2.f) * 0.5f;
      na *= persistence;
      nf *= lacunarity;
      per = scale_period(per, lacunarity);
    }

    output[index] = n;
    return;
  }

  if (noise_id == 1)
  {
    output[index] = base_perlin_fbm(pos,
                                    octaves,
                                    new_weight,
                                    persistence,
                                    lacunarity,
                                    seed,
                                    period);
  }
  else if (noise_id == 2)
  {
    output[index] = base_perlin_billow_fbm(pos,
                                           octaves,
                                           new_weight,
                                           persistence,
                                           lacunarity,
                                           seed,
                                           period);
  }
  else if (noise_id == 3)
  {
    output[index] = base_perlin_half_fbm(pos,
                                         octaves,
                                         new_weight,
                                         persistence,
                                         lacunarity,
                                         seed,
                                         period);
  }
  else if (noise_id == 4 || noise_id == 5)
  {
    output[index] = base_simplex2_fbm(pos,
                                      octaves,
                                      new_weight,
                                      persistence,
                                      lacunarity,
                                      seed);
  }
  else if (noise_id == 6)
  {
    output[index] = base_value_fbm(pos,
                                   octaves,
                                   new_weight,
                                   persistence,
                                   lacunarity,
                                   seed,
                                   period);
  }
  else if (noise_id == 7)
  {
    output[index] = base_value_cubic_fbm(pos,
                                         octaves,
                                         new_weight,
                                         persistence,
                                         lacunarity,
                                         seed,
                                         period);
  }
  else if (noise_id == 9)
  {
    output[index] = base_value_linear_fbm(pos,
                                          octaves,
                                          new_weight,
                                          persistence,
                                          lacunarity,
                                          seed,
                                          period);
  }
  else if (noise_id == 10)
  {
    output[index] = base_worley_fbm(pos,
                                    octaves,
                                    new_weight,
                                    persistence,
                                    lacunarity,
                                    seed,
                                    period);
  }
  else if (noise_id == 11)
  {
    output[index] = base_worley_double_fbm(pos,
                                           octaves,
                                           new_weight,
                                           persistence,
                                           lacunarity,
                                           seed,
                                           period);
  }
  else if (noise_id == 12)
  {
    output[index] = base_worley_value_fbm(pos,
                                          octaves,
                                          new_weight,
                                          persistence,
                                          lacunarity,
                                          seed,
                                          period);
  }
  else
  {
    output[index] = 0.f;
  }
}

void kernel noise_slope(global float *output,
                        global float *ctrl_param,
                        global float *noise_x,
                        global float *noise_y,
                        const int     nx,
                        const int     ny,
                        const int     noise_id,
                        float         kx,
                        float         ky,
                        const uint    seed,
                        const int     octaves,
                        const float   weight,
                        const float   persistence,
                        const float   lacunarity,
                        const float   gradient_scale,
                        const float   slope_power,
                        const float   octave_blend,
                        const float   bias_angle,
                        const float   bias_strength,
                        const int     has_ctrl_param,
                        const int     has_noise_x,
                        const int     has_noise_y,
                        const int     period_x,
                        const int     period_y,
                        const float4  bbox,
                        const float   angle,
                        const float   stretch,
                        const float   stretch_decay)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float ct = has_ctrl_param > 0 ? ctrl_param[index] : 1.f;
  float dx_in = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy_in = has_noise_y > 0 ? noise_y[index] : 0.f;

  int2   period = {period_x, period_y};
  float  new_weight = (1.f - ct) + ct * weight;
  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx_in, dy_in, bbox);

  float theta = angle * 0.017453292519943295f;
  float cos_a = cos(theta);
  float sin_a = sin(theta);
  float decay_base = 1.f - clamp(stretch_decay, 0.f, 1.f);

  float bias_theta = bias_angle * 0.017453292519943295f;
  float bias_cos = cos(bias_theta);
  float bias_sin = sin(bias_theta);
  float clamped_bias = clamp(bias_strength, 0.f, 1.f);

  const float eps = 0.001f;
  float2      eps_x = {eps * kx, 0.f};
  float2      eps_y = {0.f, eps * ky};

  float sum = 0.f;
  float amp = 0.6f;
  float nf = 1.f;
  float dx_sum = 0.f;
  float dy_sum = 0.f;
  int2  per = period;

  for (int k = 0; k < octaves; k++)
  {
    float2 p_s = pos * nf;
    if (stretch != 1.f)
    {
      float stretch_k = 1.f + (stretch - 1.f) * pow(decay_base, (float)k);
      p_s = apply_stretch(p_s, cos_a, sin_a, stretch_k);
    }

    float value = eval_base_noise(noise_id, p_s, seed, per);
    float dkx = (eval_base_noise(noise_id, p_s + eps_x, seed, per) -
                 eval_base_noise(noise_id, p_s - eps_x, seed, per)) /
                (2.f * eps);
    float dky = (eval_base_noise(noise_id, p_s + eps_y, seed, per) -
                 eval_base_noise(noise_id, p_s - eps_y, seed, per)) /
                (2.f * eps);

    dx_sum += amp * dkx;
    dy_sum += amp * dky;

    float grad_norm = sqrt(dx_sum * dx_sum + dy_sum * dy_sum);
    float gs = gradient_scale * grad_norm;
    float gs_p = pow(gs, slope_power);
    float w_slope = gs_p / (1.f + gs_p);

    if (clamped_bias > 0.f && grad_norm > 1e-6f)
    {
      float cos_align = (dx_sum * bias_cos + dy_sum * bias_sin) / grad_norm;
      float w_dir = (1.f - clamped_bias) +
                    clamped_bias * 0.5f * (cos_align + 1.f);
      w_slope *= w_dir;
    }

    float alpha = 1.f - exp(-octave_blend * (float)k);
    float w_octave = (1.f - alpha) + alpha * w_slope;

    sum += value * amp * w_octave;
    amp *= (1.f - new_weight) + new_weight * min(value + 1.f, 2.f) * 0.5f;
    amp *= persistence;
    nf *= lacunarity;
    per = scale_period(per, lacunarity);
  }

  output[index] = sum;
}
)""
