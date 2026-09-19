R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

// Adapted from Rune Skovbo Johansen (runevision)
// "Fast and Gorgeous Erosion Filter" / "Advanced Terrain Erosion Filter"
// https://www.shadertoy.com/view/wXcfWn
// https://blog.runevision.com/2026/03/fast-and-gorgeous-erosion-filter.html
// License: MPL 2.0 (see THIRD_PARTY_LICENSES.md)

float helper_erosion_clamp01(float v)
{
  return clamp(v, 0.f, 1.f);
}

float helper_erosion_ease_out(float t)
{
  float v = 1.f - helper_erosion_clamp01(t);
  return 1.f - v * v;
}

float helper_erosion_pow_inv(float t, float power)
{
  return 1.f - pow_float(1.f - helper_erosion_clamp01(t), power);
}

float helper_erosion_smooth_start(float t, float smoothing)
{
  if (smoothing <= 0.f) return t;
  if (t >= smoothing) return t - 0.5f * smoothing;
  return 0.5f * t * t / smoothing;
}

// Phacelle directional noise with smooth partial normalization.
// Evaluates 4x4 cellular grid with cosine and sine wave interpolation.
float4 phacelle_directional_noise(const float2 p,
                                  const float2 norm_dir,
                                  const float  freq,
                                  const float  offset_val,
                                  const float  normalization,
                                  const float  fseed)
{
  const float tau = 6.28318530718f;
  float2      side_dir = (float2)(-norm_dir.y, norm_dir.x) * freq * tau;
  float       offset = offset_val * tau;

  float2 ip = floor(p);
  float2 fp = p - ip;

  float2 phase_dir = (float2)(0.f, 0.f);
  float  weight_sum = 0.f;

  for (int i = -1; i <= 2; ++i)
  {
    for (int j = -1; j <= 2; ++j)
    {
      float2 grid_offset = (float2)((float)i, (float)j);
      float2 grid_point = ip + grid_offset;

      // Random cell offset in [-0.5, 0.5]
      float2 rnd = (hash22f_poly(grid_point, fseed) - 0.5f);
      float2 vec_from_cell = fp - grid_offset - rnd;

      float sqr_dist = dot(vec_from_cell, vec_from_cell);
      float weight = exp(-sqr_dist * 2.f);

      // Subtract to ensure weight hits zero smoothly at dist 1.5
      weight = max(0.f, weight - 0.01111f);
      weight_sum += weight;

      float wave_input = dot(vec_from_cell, side_dir) + offset;
      phase_dir += (float2)(cos(wave_input), sin(wave_input)) * weight;
    }
  }

  if (weight_sum > 1e-7f)
  {
    phase_dir /= weight_sum;
  }

  float magnitude = length(phase_dir);
  float norm_threshold = max(1.f - normalization, 1e-7f);
  magnitude = max(norm_threshold, magnitude);

  float2 norm_phase = phase_dir / magnitude;
  return (float4)(norm_phase.x, norm_phase.y, side_dir.x, side_dir.y);
}

void kernel erosion_filter(read_only image2d_t height_in,
                           global float       *output_height,
                           global float       *output_ridge_map,
                           global float       *fade_target_in,
                           const int           nx,
                           const int           ny,
                           const uint          seed,
                           const float         scale,
                           const float         strength,
                           const float         gully_weight,
                           const float         detail,
                           const float4        rounding,
                           const float4        onset,
                           const float2        assumed_slope,
                           const float         cell_scale,
                           const int           octaves,
                           const float         gain,
                           const float         lacunarity,
                           const float         normalization,
                           const float         curvature_scaling,
                           const int           has_fade_target,
                           const int           has_ridge_output,
                           const float4        bbox)
{
  int2 g = (int2)(get_global_id(0), get_global_id(1));
  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  uint  rng_state = wang_hash(seed);
  float fseed = rand(&rng_state);

  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, 0.f, 0.f, bbox);

  const sampler_t sampler_itp = CLK_NORMALIZED_COORDS_TRUE |
                                CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

  float base_h = read_imagef(height_in, sampler_itp, pos).x;

  // sample gradients and curvature from input heightmap
  float eps = 0.5f / (float)max(nx, ny);
  float h_px = read_imagef(height_in, sampler_itp, pos + (float2)(eps, 0.f)).x;
  float h_mx = read_imagef(height_in, sampler_itp, pos - (float2)(eps, 0.f)).x;
  float h_py = read_imagef(height_in, sampler_itp, pos + (float2)(0.f, eps)).x;
  float h_my = read_imagef(height_in, sampler_itp, pos - (float2)(0.f, eps)).x;

  float2 slope = (float2)((h_px - h_mx) / (2.f * eps),
                          (h_py - h_my) / (2.f * eps));

  // compute local curvature from second derivatives
  float h_p2x =
      read_imagef(height_in, sampler_itp, pos + (float2)(2.f * eps, 0.f)).x;
  float h_m2x =
      read_imagef(height_in, sampler_itp, pos - (float2)(2.f * eps, 0.f)).x;
  float h_p2y =
      read_imagef(height_in, sampler_itp, pos + (float2)(0.f, 2.f * eps)).x;
  float h_m2y =
      read_imagef(height_in, sampler_itp, pos - (float2)(0.f, 2.f * eps)).x;

  float d2x = (h_p2x + h_m2x - 2.f * base_h) / (4.f * eps * eps);
  float d2y = (h_p2y + h_m2y - 2.f * base_h) / (4.f * eps * eps);
  float curv = length((float2)(d2x, d2y)) * scale;

  // Determine fade target: [-1, 1] from valleys to peaks
  float fade_target;
  if (has_fade_target > 0)
  {
    fade_target = clamp(fade_target_in[index], -1.f, 1.f);
  }
  else
  {
    fade_target = clamp(base_h * 2.f - 1.f, -1.f, 1.f);
  }

  float3 height_and_slope = (float3)(base_h, slope.x, slope.y);
  float  cur_strength = strength * scale;
  float  freq = 1.f / (scale * cell_scale);
  float  slope_len =
      max(length(height_and_slope.yz) + curvature_scaling * curv, 1e-7f);
  float  rounding_mult = 1.f;

  float rounding_for_input = mix(rounding.y,
                                 rounding.x,
                                 helper_erosion_clamp01(fade_target + 0.5f)) *
                             rounding.z;

  float combi_mask = helper_erosion_ease_out(
      helper_erosion_smooth_start(slope_len * onset.x,
                                  rounding_for_input * onset.x));

  float ridge_combi_mask = helper_erosion_ease_out(slope_len * onset.z);
  float ridge_fade_target = fade_target;

  // Mix actual slope with assumed slope for initial gully direction
  float2 gully_slope = mix(height_and_slope.yz,
                           height_and_slope.yz / slope_len * assumed_slope.x,
                           assumed_slope.y);

  for (int oct = 0; oct < octaves; ++oct)
  {
    float2 norm_gully_dir = normed2(gully_slope);
    if (dot(norm_gully_dir, norm_gully_dir) < 1e-7f)
    {
      norm_gully_dir = (float2)(0.f, 1.f);
    }

    float4 phacelle = phacelle_directional_noise(pos * freq,
                                                 norm_gully_dir,
                                                 cell_scale,
                                                 0.25f,
                                                 normalization,
                                                 fseed + (float)oct * 0.13f);

    phacelle.zw *= -freq;
    float sloping = fabs(phacelle.y);

    // Straight-gullies slope accumulation
    float sign_val = phacelle.y < 0.f ? -1.f : 1.f;
    gully_slope += sign_val * phacelle.zw * cur_strength * gully_weight;

    float3 gullies = (float3)(phacelle.x,
                              phacelle.y * phacelle.z,
                              phacelle.y * phacelle.w);
    float3 faded_gullies = mix((float3)(fade_target, 0.f, 0.f),
                               gullies * gully_weight,
                               combi_mask);

    height_and_slope += faded_gullies * cur_strength;
    fade_target = faded_gullies.x;

    float rounding_for_octave = mix(rounding.y,
                                    rounding.x,
                                    helper_erosion_clamp01(phacelle.x + 0.5f)) *
                                rounding_mult;

    float new_mask = helper_erosion_ease_out(
        helper_erosion_smooth_start(sloping * onset.y,
                                    rounding_for_octave * onset.y));
    combi_mask = helper_erosion_pow_inv(combi_mask, detail) * new_mask;

    // Ridge map accumulation
    ridge_fade_target = mix(ridge_fade_target, gullies.x, ridge_combi_mask);
    float new_ridge_mask = helper_erosion_ease_out(sloping * onset.w);
    ridge_combi_mask *= new_ridge_mask;

    cur_strength *= gain;
    freq *= lacunarity;
    rounding_mult *= rounding.w;
  }

  output_height[index] = height_and_slope.x;

  if (has_ridge_output > 0)
  {
    output_ridge_map[index] = ridge_fade_target * (1.f - ridge_combi_mask);
  }
}
)""
