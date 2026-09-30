R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel voronoi(global float *output,
                    global float *ctrl_param,
                    global float *noise_x,
                    global float *noise_y,
                    const int     nx,
                    const int     ny,
                    const float   kx,
                    const float   ky,
                    const uint    seed,
                    const float2  jitter,
                    const float   k_smoothing,
                    const float   exp_sigma,
                    const int     return_type,
                    const int     has_ctrl_param,
                    const int     has_noise_x,
                    const int     has_noise_y,
                    const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float ct = has_ctrl_param > 0 ? ctrl_param[index] : 1.f;
  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  float val = 0.f;

  switch (return_type)
  {
  case 0: val = base_voronoi_f1(pos, ct * jitter, k_smoothing, fseed); break;
  case 1: val = base_voronoi_f2(pos, ct * jitter, k_smoothing, fseed); break;
  case 2: val = base_voronoi_f1tf2(pos, ct * jitter, k_smoothing, fseed); break;
  case 3: val = base_voronoi_f1df2(pos, ct * jitter, k_smoothing, fseed); break;
  case 4: val = base_voronoi_f2mf1(pos, ct * jitter, k_smoothing, fseed); break;
  case 5:
    val = base_voronoi_edge_distance_exp(pos,
                                         ct * jitter,
                                         k_smoothing,
                                         exp_sigma,
                                         fseed);
    break;
  case 6:
    val = base_voronoi_edge_distance(pos, ct * jitter, k_smoothing, fseed);
    break;
  case 7:
    val = base_voronoi_constant(pos, ct * jitter, k_smoothing, fseed);
    break;
  case 8:
    val = base_voronoi_constant_f2mf1(pos, ct * jitter, k_smoothing, fseed);
    break;
  }

  output[index] = val;
}

void kernel voronoi_fbm(global float *output,
                        global float *ctrl_param,
                        global float *noise_x,
                        global float *noise_y,
                        const int     nx,
                        const int     ny,
                        const float   kx,
                        const float   ky,
                        const uint    seed,
                        const float2  jitter,
                        const float   k_smoothing,
                        const float   exp_sigma,
                        const int     return_type,
                        const int     octaves,
                        const float   weight,
                        const float   persistence,
                        const float   lacunarity,
                        const int     has_ctrl_param,
                        const int     has_noise_x,
                        const int     has_noise_y,
                        const float4  bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float ct = has_ctrl_param > 0 ? ctrl_param[index] : 1.f;
  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, kx, ky, dx, dy, bbox);

  float val = 0.f;

  switch (return_type)
  {
  case 0:
    val = base_voronoi_f1_fbm(pos,
                              ct * jitter,
                              k_smoothing,
                              octaves,
                              weight,
                              persistence,
                              lacunarity,
                              fseed);
    break;
  case 1:
    val = base_voronoi_f2_fbm(pos,
                              ct * jitter,
                              k_smoothing,
                              octaves,
                              weight,
                              persistence,
                              lacunarity,
                              fseed);
    break;
  case 2:
    val = base_voronoi_f1tf2_fbm(pos,
                                 ct * jitter,
                                 k_smoothing,
                                 octaves,
                                 weight,
                                 persistence,
                                 lacunarity,
                                 fseed);
    break;
  case 3:
    val = base_voronoi_f1df2_fbm(pos,
                                 ct * jitter,
                                 k_smoothing,
                                 octaves,
                                 weight,
                                 persistence,
                                 lacunarity,
                                 fseed);
    break;
  case 4:
    val = base_voronoi_f2mf1_fbm(pos,
                                 ct * jitter,
                                 k_smoothing,
                                 octaves,
                                 weight,
                                 persistence,
                                 lacunarity,
                                 fseed);
    break;
  case 5:
    val = base_voronoi_edge_distance_exp_fbm(pos,
                                             ct * jitter,
                                             k_smoothing,
                                             exp_sigma,
                                             octaves,
                                             weight,
                                             persistence,
                                             lacunarity,
                                             fseed);
    break;
  case 6:
    val = base_voronoi_edge_distance_fbm(pos,
                                         ct * jitter,
                                         k_smoothing,
                                         octaves,
                                         weight,
                                         persistence,
                                         lacunarity,
                                         fseed);
    break;
  case 7:
    val = base_voronoi_constant_fbm(pos,
                                    ct * jitter,
                                    k_smoothing,
                                    octaves,
                                    weight,
                                    persistence,
                                    lacunarity,
                                    fseed);
    break;
  case 8:
    val = base_voronoi_constant_f2mf1_fbm(pos,
                                          ct * jitter,
                                          k_smoothing,
                                          octaves,
                                          weight,
                                          persistence,
                                          lacunarity,
                                          fseed);
    break;
  }

  output[index] = val;
}
)""
