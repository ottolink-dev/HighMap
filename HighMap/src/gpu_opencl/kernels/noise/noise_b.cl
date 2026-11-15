R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// --- OCL KERNEL

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
                  const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  int2 period = {period_x, period_y};

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  if (noise_id == 1)
  {
    output[index] = base_perlin(pos, seed, period);
  }
  else if (noise_id == 2)
  {
    float value = base_perlin(pos, seed, period);
    output[index] = 2.f * fabs(value) - 1.f;
  }
  else if (noise_id == 3)
  {
    float value = base_perlin(pos, seed, period);
    output[index] = max_smooth(value, 0.f, 0.5f);
  }
  else if (noise_id == 4 || noise_id == 5)
  {
    output[index] = base_simplex2(pos, seed);
  }
  else if (noise_id == 6)
  {
    output[index] = base_value(pos, seed, period);
  }
  else if (noise_id == 7)
  {
    output[index] = base_value_cubic(pos, seed, period);
  }
  else if (noise_id == 9)
  {
    output[index] = base_value_linear(pos, seed, period);
  }
  else if (noise_id == 10)
  {
    output[index] = base_worley(pos, seed, period);
  }
  else if (noise_id == 11)
  {
    output[index] = base_worley_double(pos, seed, 0.f, 2.f, period);
  }
  else if (noise_id == 12)
  {
    output[index] = base_worley_value(pos, seed, period);
  }
  else
  {
    output[index] = 0.f;
  }
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
                      const float4  bbox)
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
)""
