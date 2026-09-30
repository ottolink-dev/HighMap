R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel rifts(global float *output,
                  global float *noise_x,
                  global float *noise_y,
                  global float *mask,
                  const int     nx,
                  const int     ny,
                  const float2  kw,
                  const float   angle,
                  const float   amplitude,
                  const uint    seed,
                  const float   elevation_noise_shift,
                  const float   k_smooth_bottom,
                  const float   k_smooth_top,
                  const float   radial_spread_amp,
                  const float   elevation_noise_amp,
                  const float   clamp_vmin,
                  const float   remap_vmin,
                  const int     apply_mask,
                  const int     reverse_mask,
                  const float   mask_gamma,
                  const int     has_noise_x,
                  const int     has_noise_y,
                  const int     has_mask,
                  const float2  center,
                  const float4  bbox)
{

  // ---

  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float dy = has_noise_y > 0 ? noise_y[index] : 0.f;

  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, dx, dy, bbox);

  // parameters
  const float val0 = output[index];
  float       val = val0;
  float       alpha = angle / 180.f * 3.141592f;
  float2      dir = (float2)(cos(alpha), sin(alpha));

  // rotate
  float dr_x = dot(pos - center.x, (float2)(cos(alpha), sin(alpha)));
  float dr_y = dot(pos - center.y, (float2)(-sin(alpha), cos(alpha)));

  float dist_to_center = length((float2)(dr_x, dr_y));
  // use a sigmoid as a smooth 'sign' operator
  dist_to_center *= -2.f / (1.f + exp(-20.f * dr_x)) + 1.f;

  float2 pos_r = (float2)(kw.x * (dr_x + radial_spread_amp * dist_to_center +
                                  elevation_noise_amp *
                                      (val0 - elevation_noise_shift)),
                          kw.y * dr_y);

  float noise_v = base_voronoi_f1df2(pos_r,
                                     (float2)(1.f, 1.f),
                                     k_smooth_bottom,
                                     fseed);
  noise_v = smax(min(noise_v, 1.f), clamp_vmin, k_smooth_top);
  noise_v = remap_from(noise_v, 1.f, remap_vmin, clamp_vmin, 1.f);

  val += amplitude * noise_v;

  // elevation mask
  {
    float t = apply_mask ? pow(val0, mask_gamma) : 1.f;
    if (reverse_mask) t = 1.f - t;
    val = lerp(val0, val, t);
  }

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
