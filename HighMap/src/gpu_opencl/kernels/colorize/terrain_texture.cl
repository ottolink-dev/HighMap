R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

float3 terrain_shade_color(float3 color, float t)
{
  t = clamp(t, 0.f, 1.f);
  const float  luma = dot(color, (float3)(0.299f, 0.587f, 0.114f));
  const float3 grey = (float3)(luma, luma, luma);
  const float3 lit = color * (1.f + 0.4f * t);
  return mix(grey * (1.f + 0.4f * t), lit, 1.f - 0.3f * t);
}

void kernel terrain_fallline_walk(read_only image2d_t  dir_x,
                                  read_only image2d_t  dir_y,
                                  read_only image2d_t  in_source,
                                  read_only image2d_t  in_n_rill,
                                  write_only image2d_t out_rill,
                                  write_only image2d_t out_streak,
                                  const int            nx,
                                  const int            ny,
                                  const int            rill_steps,
                                  const int            streak_steps)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler_nn = CLK_NORMALIZED_COORDS_FALSE |
                               CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float rill_sum = read_imagef(in_n_rill, sampler_nn, g).x;
  int   rill_count = 1;

  float streak_best = read_imagef(in_source, sampler_nn, g).x;

  const int max_steps = max(rill_steps, streak_steps);

  float2 cur = (float2)((float)g.x, (float)g.y);

  for (int s = 1; s <= max_steps; ++s)
  {
    int2 coord = (int2)((int)(cur.x + 0.5f), (int)(cur.y + 0.5f));

    float gx = read_imagef(dir_x, sampler_nn, coord).x;
    float gy = read_imagef(dir_y, sampler_nn, coord).x;

    if (fabs(gx) < 1e-6f && fabs(gy) < 1e-6f)
    {
      break;
    }

    cur.x += gx;
    cur.y += gy;

    if (cur.x < 0.f || cur.y < 0.f || cur.x > (float)(nx - 1) ||
        cur.y > (float)(ny - 1))
    {
      break;
    }

    int2 next_coord = (int2)((int)(cur.x + 0.5f), (int)(cur.y + 0.5f));

    if (s <= rill_steps)
    {
      rill_sum += read_imagef(in_n_rill, sampler_nn, next_coord).x;
      rill_count++;
    }

    if (s <= streak_steps)
    {
      float src_val = read_imagef(in_source, sampler_nn, next_coord).x;
      float decay = 1.f - (float)s / (float)(streak_steps + 1);
      streak_best = max(streak_best, src_val * decay);
    }
  }

  float rill_res = clamp(rill_sum / sqrt((float)rill_count), -1.f, 1.f);

  write_imagef(out_rill, g, rill_res);
  write_imagef(out_streak, g, streak_best);
}

void kernel terrain_shade(read_only image2d_t in_z,
                          read_only image2d_t in_slope,
                          read_only image2d_t in_crest,
                          read_only image2d_t in_region,
                          read_only image2d_t in_fine,
                          read_only image2d_t in_flow,
                          read_only image2d_t in_deposit,
                          read_only image2d_t in_moisture,
                          read_only image2d_t in_streak,
                          read_only image2d_t in_rill,
                          read_only image2d_t in_nl,
                          read_only image2d_t in_nm,
                          read_only image2d_t in_nf,
                          read_only image2d_t in_ng,
                          global float       *out_r,
                          global float       *out_g,
                          global float       *out_b,
                          global float       *out_a,
                          const int           nx,
                          const int           ny,
                          const float         zmin,
                          const float         zrange,
                          const float         rock_talus,
                          const float         rock_edge,
                          const float         strata,
                          const float         snow_line,
                          const float         snow_amount,
                          const float         color_noise,
                          const float         flow_width,
                          const float         param_flow,
                          const float         param_deposition,
                          const float         param_streaks,
                          const float         param_cavity,
                          const float         param_detail,
                          const float3        col_grass,
                          const float3        col_dry,
                          const float3        col_soil,
                          const float3        col_sediment,
                          const float3        col_rock,
                          const float3        col_ridge,
                          const float3        col_snow)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const int idx = g.y * nx + g.x;

  const sampler_t sampler_nn = CLK_NORMALIZED_COORDS_FALSE |
                               CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const float zv = read_imagef(in_z, sampler_nn, g).x;
  const float h = (zv - zmin) / zrange;
  const float slope_val = read_imagef(in_slope, sampler_nn, g).x;
  const float s = slope_val / rock_talus;
  const float c = read_imagef(in_crest, sampler_nn, g).x;
  const float r = read_imagef(in_region, sampler_nn, g).x;
  const float fine_val = read_imagef(in_fine, sampler_nn, g).x;
  const float flow_val = read_imagef(in_flow, sampler_nn, g).x;
  const float dep_val = read_imagef(in_deposit, sampler_nn, g).x;
  const float moist_val = read_imagef(in_moisture, sampler_nn, g).x;
  const float streak_val_raw = read_imagef(in_streak, sampler_nn, g).x;
  const float rl = read_imagef(in_rill, sampler_nn, g).x;

  const float nl = read_imagef(in_nl, sampler_nn, g).x;
  const float nm = read_imagef(in_nm, sampler_nn, g).x;
  const float nf = read_imagef(in_nf, sampler_nn, g).x;
  const float ng = read_imagef(in_ng, sampler_nn, g).x;

  const float cn = color_noise;
  const float two_pi = 6.283185307179586f;

  // --- Grass & dry grass
  const float wet_val = 0.5f + 0.6f * moist_val - 0.35f * r - 0.2f * c +
                        0.15f * nl;
  const float dry_val = 1.f - wet_val + 0.35f * (h - 0.4f) + 0.12f * nm +
                        0.08f * nf;
  const float dry_mask = smoothstep3_gl(dry_val, 0.42f, 0.58f);

  const float  grass_sh = cn * (0.3f * nm + 0.3f * nf + 0.2f * ng);
  const float3 c_grass = terrain_shade_color(col_grass, grass_sh);
  const float3 c_dry = terrain_shade_color(col_dry, grass_sh);

  float3 col = mix(c_grass, c_dry, dry_mask);

  // --- Soil
  const float soil_val = s + 0.2f * nm + 0.1f * nf + 0.25f * max(c, 0.f);
  const float soil_mask = smoothstep3_gl(soil_val, 0.35f, 0.6f);
  const float soil_sh = cn * (0.25f * rl + 0.3f * nm + 0.15f * ng);

  col = mix(col, terrain_shade_color(col_soil, soil_sh), soil_mask);

  // --- Sediment
  const float flow_mask = smoothstep3_gl(flow_val + 0.05f * nf,
                                         flow_width,
                                         flow_width + 0.08f);
  const float dep_mask = smoothstep3_gl(dep_val + 0.15f * nf, 0.4f, 0.8f);

  float sed_mask = max(param_flow * flow_mask, param_deposition * dep_mask);
  sed_mask *= 1.f - smoothstep3_gl(s, 0.6f, 1.f);

  const float sed_sh = cn * (0.4f * nf + 0.2f * nm + 0.2f * ng);
  col = mix(col, terrain_shade_color(col_sediment, sed_sh), sed_mask);

  // --- Rock & strata
  const float rock_val = s + 0.25f * nm + 0.1f * nf + 0.3f * max(c, 0.f) * h;
  const float rock_mask = smoothstep3_gl(rock_val,
                                         1.f - rock_edge,
                                         1.f + rock_edge);

  const float band = sin(two_pi * (14.f * h + 0.4f * nl + 0.1f * nm));
  const float rock_sh = 0.35f * strata * band +
                        cn * (0.2f * rl + 0.3f * nf + 0.15f * ng) + 0.25f * c;

  col = mix(col, terrain_shade_color(col_rock, rock_sh), rock_mask);

  // --- Ridge streaks
  const float streak_val = streak_val_raw * (0.75f + 0.25f * (1.f + nf));
  const float streak_mask = smoothstep3_gl(streak_val, 0.25f, 0.65f);
  const float slope_gate = smoothstep3_gl(s, 0.15f, 0.4f);
  const float spill_mask = param_streaks * streak_mask * slope_gate;

  const float ridge_sh = cn * (0.3f * nf + 0.15f * rl + 0.15f * ng);
  col = mix(col, terrain_shade_color(col_ridge, ridge_sh), spill_mask);

  // --- Snow
  const float snow_val = h - snow_line + 0.08f * nm + 0.04f * nf -
                         0.15f * max(0.f, s - 0.9f) - 0.05f * c;
  const float snow_mask = snow_amount * smoothstep3_gl(snow_val, -0.04f, 0.04f);
  const float snow_sh = 0.15f * nf - 0.3f * max(-c, 0.f);

  col = mix(col, terrain_shade_color(col_snow, snow_sh), snow_mask);

  // --- Micro-relief & ambient shading
  const float cavity_factor = clamp(-0.7f * c - 0.4f * r, 0.f, 1.f);
  col *= 1.f - param_cavity * 0.35f * cavity_factor;
  col *= 1.f + param_detail * 0.22f * fine_val;
  col *= 1.f + cn * (0.05f * nm + 0.06f * nf + 0.05f * ng);
  col += cn * 0.015f * nl * (float3)(1.f, 0.6f, 0.1f);

  // --- Output RGBA planar arrays
  out_r[idx] = clamp(col.x, 0.f, 1.f);
  out_g[idx] = clamp(col.y, 0.f, 1.f);
  out_b[idx] = clamp(col.z, 0.f, 1.f);
  out_a[idx] = 1.f;
}
)""
